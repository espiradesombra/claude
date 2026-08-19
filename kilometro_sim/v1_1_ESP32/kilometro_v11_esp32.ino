/*
 * Kilometro prototipo v1.1 — ESP32
 * Bateria de reset por lastre de perneado
 *
 * - Pernos: SOL_O, SOL_RA (R alta), SOL_RB (R baja)
 * - Make-before-break en enganche y entrega
 * - FC_ALTA / FC_BAJA interlocks
 * - Log CSV por Serial 115200
 *
 * Arduino IDE: placa "ESP32 Dev Module"
 * No alimentar solenoides desde GPIO.
 */

#include <Arduino.h>

// -------------------- Pines --------------------
static const int PIN_SOL_O  = 25;
static const int PIN_SOL_RA = 26;
static const int PIN_SOL_RB = 27;
static const int PIN_FC_ALTA = 32;  // INPUT_PULLUP, activo LOW
static const int PIN_FC_BAJA = 33;
static const int PIN_BTN_START = 18;  // pull-up interno
static const int PIN_BTN_RESET = 19;
static const int PIN_LED_OK    = 2;
static const int PIN_LED_FAULT = 4;

// true = activo en LOW (microswitch a GND)
static const bool FC_ACTIVE_LOW = true;
static const bool BTN_ACTIVE_LOW = true;

// -------------------- Tiempos (ms) --------------------
static const uint32_t HOLD_MS = 300;
static const uint32_t SETTLE_MS = 150;
static const uint32_t SOL_PULSE_MS = 250;      // pulso biestable
static const uint32_t T_BAJADA_MAX_MS = 30000;
static const uint32_t T_SUBIDA_MAX_MS = 45000;
static const uint32_t DEBOUNCE_MS = 40;

// Si solenoides necesitan quedar ON (no biestables), poner true
// y respetar duty del fabricante.
static const bool SOLENOID_HOLD_MODE = false;

// -------------------- FSM --------------------
enum State : uint8_t {
  ST_IDLE_ALTA = 0,
  ST_ENGANCHE,
  ST_ESPERA_BAJADA,
  ST_ENTREGA,
  ST_ESPERA_SUBIDA,
  ST_CICLO_OK,
  ST_FAULT
};

static State state = ST_IDLE_ALTA;
static uint32_t stateEnterMs = 0;
static uint32_t tMark = 0;

static uint32_t cycleCount = 0;
static uint32_t tEnganche = 0, tBajada = 0, tEntrega = 0, tSubida = 0, tTotal = 0;
static uint32_t cycleStartMs = 0;
static int faultCode = 0;

// Software state of latches (1 = locked / engaged)
static bool latO = false;
static bool latRA = true;   // en idle, lastre en guia alta por defecto
static bool latRB = false;

// -------------------- Utilidades IO --------------------
static bool readActive(int pin, bool activeLow) {
  int v = digitalRead(pin);
  return activeLow ? (v == LOW) : (v == HIGH);
}

static bool fcAlta() { return readActive(PIN_FC_ALTA, FC_ACTIVE_LOW); }
static bool fcBaja() { return readActive(PIN_FC_BAJA, FC_ACTIVE_LOW); }

static bool btnStartPressed() {
  static uint32_t last = 0;
  static bool lastState = false;
  bool now = readActive(PIN_BTN_START, BTN_ACTIVE_LOW);
  if (now && !lastState && (millis() - last > DEBOUNCE_MS)) {
    last = millis();
    lastState = now;
    return true;
  }
  lastState = now;
  return false;
}

static bool btnResetPressed() {
  static uint32_t last = 0;
  static bool lastState = false;
  bool now = readActive(PIN_BTN_RESET, BTN_ACTIVE_LOW);
  if (now && !lastState && (millis() - last > DEBOUNCE_MS)) {
    last = millis();
    lastState = now;
    return true;
  }
  lastState = now;
  return false;
}

static void solWrite(int pin, bool on) {
  digitalWrite(pin, on ? HIGH : LOW);
}

// Pulso o hold segun modo
static void solPulse(int pin, bool /*engage*/) {
  solWrite(pin, true);
  delay(SOL_PULSE_MS);
  if (!SOLENOID_HOLD_MODE) {
    solWrite(pin, false);
  }
}

static void enter(State s) {
  state = s;
  stateEnterMs = millis();
}

static void setFault(int code) {
  faultCode = code;
  enter(ST_FAULT);
  digitalWrite(PIN_LED_FAULT, HIGH);
  digitalWrite(PIN_LED_OK, LOW);
  // No abrir ambos lados: apagar drives (pasadores biestables mantienen muesca)
  solWrite(PIN_SOL_O, false);
  solWrite(PIN_SOL_RA, false);
  solWrite(PIN_SOL_RB, false);
  Serial.print(F("FAULT code="));
  Serial.println(code);
}

static void printHelp() {
  Serial.println(F("Kilometro v1.1 ESP32"));
  Serial.println(F("s = start cycle"));
  Serial.println(F("r = reset fault"));
  Serial.println(F("p = print state"));
  Serial.println(F("z = zero cycle counter"));
  Serial.println(F("h = help"));
  Serial.println(F("CSV: cycle,t_enganche_ms,t_bajada_ms,t_entrega_ms,t_subida_ms,t_total_ms,ok,fault_code"));
}

static const char* stateName(State s) {
  switch (s) {
    case ST_IDLE_ALTA: return "IDLE_ALTA";
    case ST_ENGANCHE: return "ENGANCHE";
    case ST_ESPERA_BAJADA: return "ESPERA_BAJADA";
    case ST_ENTREGA: return "ENTREGA";
    case ST_ESPERA_SUBIDA: return "ESPERA_SUBIDA";
    case ST_CICLO_OK: return "CICLO_OK";
    case ST_FAULT: return "FAULT";
    default: return "?";
  }
}

static void printState() {
  Serial.print(F("state="));
  Serial.print(stateName(state));
  Serial.print(F(" fcA="));
  Serial.print(fcAlta() ? 1 : 0);
  Serial.print(F(" fcB="));
  Serial.print(fcBaja() ? 1 : 0);
  Serial.print(F(" latO="));
  Serial.print(latO ? 1 : 0);
  Serial.print(F(" latRA="));
  Serial.print(latRA ? 1 : 0);
  Serial.print(F(" latRB="));
  Serial.print(latRB ? 1 : 0);
  Serial.print(F(" cycles="));
  Serial.print(cycleCount);
  Serial.print(F(" fault="));
  Serial.println(faultCode);
}

static void logCsv(bool ok) {
  Serial.print(cycleCount);
  Serial.print(',');
  Serial.print(tEnganche);
  Serial.print(',');
  Serial.print(tBajada);
  Serial.print(',');
  Serial.print(tEntrega);
  Serial.print(',');
  Serial.print(tSubida);
  Serial.print(',');
  Serial.print(tTotal);
  Serial.print(',');
  Serial.print(ok ? 1 : 0);
  Serial.print(',');
  Serial.println(faultCode);
}

// -------------------- Secuencias pernos --------------------
// Enganche: lastre de guia alta -> objeto
// Make-before-break: O ON, luego RA OFF
static bool doEnganche() {
  if (!fcAlta()) {
    setFault(1);  // no en ALTA
    return false;
  }
  uint32_t t0 = millis();

  // Cerrar O (sujeta lastre al modulo)
  solPulse(PIN_SOL_O, true);
  latO = true;
  delay(HOLD_MS);

  // Abrir R_alta (libera de guia)
  solPulse(PIN_SOL_RA, true);  // pulso en sentido "open" si es biestable
  // En hold mode: energizar RA off coil — aqui un solo solenoide "toggle"
  if (SOLENOID_HOLD_MODE) {
    solWrite(PIN_SOL_RA, false);  // si RA ON = locked, apagar = soltar (depende montaje)
  }
  latRA = false;
  delay(SETTLE_MS);

  tEnganche = millis() - t0;
  return true;
}

// Entrega: lastre de objeto -> guia baja
// Make-before-break: RB ON, luego O OFF
static bool doEntrega() {
  if (!fcBaja()) {
    setFault(2);  // no en BAJA
    return false;
  }
  uint32_t t0 = millis();

  solPulse(PIN_SOL_RB, true);
  latRB = true;
  delay(HOLD_MS);

  solPulse(PIN_SOL_O, true);  // pulso a posicion open si biestable
  if (SOLENOID_HOLD_MODE) {
    solWrite(PIN_SOL_O, false);
  }
  latO = false;
  delay(SETTLE_MS);

  tEntrega = millis() - t0;
  return true;
}

// -------------------- Setup / Loop --------------------
void setup() {
  Serial.begin(115200);
  delay(200);

  pinMode(PIN_SOL_O, OUTPUT);
  pinMode(PIN_SOL_RA, OUTPUT);
  pinMode(PIN_SOL_RB, OUTPUT);
  solWrite(PIN_SOL_O, false);
  solWrite(PIN_SOL_RA, false);
  solWrite(PIN_SOL_RB, false);

  pinMode(PIN_FC_ALTA, INPUT_PULLUP);
  pinMode(PIN_FC_BAJA, INPUT_PULLUP);
  pinMode(PIN_BTN_START, INPUT_PULLUP);
  pinMode(PIN_BTN_RESET, INPUT_PULLUP);
  pinMode(PIN_LED_OK, OUTPUT);
  pinMode(PIN_LED_FAULT, OUTPUT);
  digitalWrite(PIN_LED_OK, HIGH);
  digitalWrite(PIN_LED_FAULT, LOW);

  printHelp();
  printState();
  Serial.println(F("Ready. Place module at ALTA, ballast on R_alta, send 's'"));
}

void loop() {
  // Serial
  if (Serial.available()) {
    char c = Serial.read();
    if (c == 'h') printHelp();
    else if (c == 'p') printState();
    else if (c == 'z') {
      cycleCount = 0;
      Serial.println(F("cycle counter zeroed"));
    } else if (c == 'r') {
      if (state == ST_FAULT) {
        faultCode = 0;
        digitalWrite(PIN_LED_FAULT, LOW);
        digitalWrite(PIN_LED_OK, HIGH);
        enter(ST_IDLE_ALTA);
        Serial.println(F("fault cleared -> IDLE_ALTA (verify mechanics!)"));
      }
    } else if (c == 's') {
      if (state == ST_IDLE_ALTA) {
        cycleStartMs = millis();
        tEnganche = tBajada = tEntrega = tSubida = tTotal = 0;
        enter(ST_ENGANCHE);
        Serial.println(F("START cycle"));
      } else {
        Serial.println(F("busy or fault"));
      }
    }
  }

  if (btnResetPressed() && state == ST_FAULT) {
    faultCode = 0;
    digitalWrite(PIN_LED_FAULT, LOW);
    digitalWrite(PIN_LED_OK, HIGH);
    enter(ST_IDLE_ALTA);
  }

  if (btnStartPressed() && state == ST_IDLE_ALTA) {
    cycleStartMs = millis();
    tEnganche = tBajada = tEntrega = tSubida = tTotal = 0;
    enter(ST_ENGANCHE);
  }

  // FSM
  switch (state) {
    case ST_IDLE_ALTA:
      digitalWrite(PIN_LED_OK, (millis() / 500) % 2);
      break;

    case ST_ENGANCHE:
      if (!doEnganche()) break;
      tMark = millis();
      enter(ST_ESPERA_BAJADA);
      Serial.println(F("enganchado -> espera bajada (FC_BAJA)"));
      break;

    case ST_ESPERA_BAJADA:
      if (fcBaja()) {
        tBajada = millis() - tMark;
        enter(ST_ENTREGA);
      } else if (millis() - stateEnterMs > T_BAJADA_MAX_MS) {
        setFault(3);  // timeout bajada
      }
      break;

    case ST_ENTREGA:
      if (!doEntrega()) break;
      tMark = millis();
      enter(ST_ESPERA_SUBIDA);
      Serial.println(F("entregado -> espera subida (FC_ALTA)"));
      break;

    case ST_ESPERA_SUBIDA:
      if (fcAlta()) {
        tSubida = millis() - tMark;
        enter(ST_CICLO_OK);
      } else if (millis() - stateEnterMs > T_SUBIDA_MAX_MS) {
        setFault(4);  // timeout subida
      }
      break;

    case ST_CICLO_OK:
      cycleCount++;
      tTotal = millis() - cycleStartMs;
      // Tras entrega, lastre en RB; para el siguiente ciclo el operador
      // debe devolver lastre a RA (recarga bateria) y resetear latches logicos:
      latRB = true;
      latRA = false;
      latO = false;
      logCsv(true);
      Serial.println(F("CICLO_OK — devuelve lastre a ALTA manualmente, luego 's'"));
      digitalWrite(PIN_LED_OK, HIGH);
      enter(ST_IDLE_ALTA);
      break;

    case ST_FAULT:
      digitalWrite(PIN_LED_FAULT, (millis() / 200) % 2);
      break;
  }
}
