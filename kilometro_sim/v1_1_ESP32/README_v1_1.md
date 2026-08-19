# Kilómetro prototipo v1.1 — Solenoides + ESP32 log

Extiende el MVP de tanque (`PROTOTIPO_TANQUE_PLANOS_Y_BOM.md`) con:

1. **Actuación de pernos** por solenoide (o servo/linear) en **O**, **R_alta**, **R_baja**
2. **Secuencia make-before-break** en firmware
3. **ESP32** que registra tiempos de ciclo en Serial (y opcional CSV por USB)
4. **Interlocks** básicos (fin de carrera ALTA/BAJA) antes de conmutar

**Sigue siendo batería de lastre**, no generador perpetuo.  
La recarga del lastre a cota ALTA puede seguir siendo **manual** en v1.1; el firmware solo automatiza enganche/entrega y el log.

---

## 1. Arquitectura

```text
  [Botones / Serial cmds]
           │
        ESP32
     ┌─────┴─────┐
     │  lógica   │
     │  FSM      │
     └─────┬─────┘
           │ GPIO → drivers MOSFET
     ┌─────┼─────┬──────────┐
     ▼     ▼     ▼          ▼
   SOL_O  SOL_RA SOL_RB   (opcional LED/buzzer)
     │     │     │
   perno  R_alta R_baja
     O
           │
     [FC_ALTA] [FC_BAJA]  ← finales de carrera (NC/NO)
```

### Estados FSM

| Estado | Significado |
|--------|-------------|
| `IDLE_ALTA` | Módulo arriba, n=3, listo |
| `ENGANCHE` | Make-before-break: O on → R_alta off |
| `ESPERA_BAJADA` | Lastre en módulo; espera FC_BAJA |
| `ENTREGA` | Make-before-break: R_baja on → O off |
| `ESPERA_SUBIDA` | Módulo ligero; espera FC_ALTA |
| `CICLO_OK` | Log del ciclo; vuelve a IDLE_ALTA |
| `FAULT` | Timeout o interlock; suelta a seguro |

### Make-before-break (firmware)

**Enganche ALTA:**

```text
assert FC_ALTA
SOL_O = ON
delay HOLD_MS
SOL_RA = OFF          // libera lastre de guía alta
delay SETTLE_MS
→ ESPERA_BAJADA
```

**Entrega BAJA:**

```text
assert FC_BAJA
SOL_RB = ON
delay HOLD_MS
SOL_O = OFF           // libera lastre del módulo
delay SETTLE_MS
→ ESPERA_SUBIDA
```

**Estado seguro en FAULT:**

```text
// Preferir lastre retenido por R si se conoce cota; si no, no abrir ambos
// Política v1.1: no apagar todos a la vez. Mantener último R/O conocidos.
```

---

## 2. Hardware

### 2.1 Electrónica (BOM adicional v1.1)

| # | Ítem | Spec | Cant. | Est. € |
|---|------|------|-------|--------|
| E1 | ESP32 DevKit | ESP32-WROOM | 1 | 6–12 |
| E2 | Fuente 12 V | 2–3 A, conmutada | 1 | 10–20 |
| E3 | Buck 12→5 V | para ESP32 si no USB | 1 | 3–6 |
| E4 | MOSFET logic-level | IRLZ44N / IRLZ34 o módulo 3 ch | 3 | 5–12 |
| E5 | Diodo flyback | 1N5408 o UF4007 en paralelo bobina | 3 | 1–2 |
| E6 | Solenoide push-pull | 12 V, carrera 5–10 mm, duty intermitente | 3 | 25–60 |
| E7 | Final de carrera | microswitch palanca IP54 o reed + imán | 2 | 4–10 |
| E8 | Botón START / RESET | momentáneo | 2 | 2 |
| E9 | LED estado | 5 mm + R 220 Ω | 2 | 1 |
| E10 | Protoboard + cables | — | 1 | 8–15 |
| E11 | Caja estanca electrónica | IP54, **fuera del tanque** | 1 | 10–20 |
| E12 | Prensaestopas | cables a solenoides | 3–6 | 8–15 |

**Total extra v1.1:** ~**80–170 €** sobre la v1 mecánica.

### 2.2 Notas de solenoides bajo agua

Los solenoides **no van sumergidos** si puedes evitarlo:

- Actuar **pasadores con bowden** (cable) desde fuera del tanque (recomendado), o  
- Usar carcasas IP68 + fuelle (más caro y frágil).

**Recomendación v1.1:** solenoides en el **puente seco** + cable bowden a cada pasador (igual que palancas manuales de v1, pero motorizadas).

### 2.3 Pines GPIO sugeridos (ESP32)

| Señal | GPIO | Notas |
|-------|------|-------|
| SOL_O | 25 | MOSFET low-side |
| SOL_RA | 26 | R_alta |
| SOL_RB | 27 | R_baja |
| FC_ALTA | 32 | INPUT_PULLUP, activo LOW al activar |
| FC_BAJA | 33 | INPUT_PULLUP |
| BTN_START | 34 | solo input (o 0 con pull-up externo) |
| BTN_RESET | 35 | solo input |
| LED_OK | 2 | onboard o externo |
| LED_FAULT | 4 | |

Ajustar si el DevKit usa otros pines ocupados.

### 2.4 Esquema eléctrico (lógica)

```text
12V+ ──┬── bobina SOL_x ── MOSFET drain
       │                      source ── GND
       └── diodo catodo a 12V+, anodo a drain (flyback)

ESP32 GPIO ── 220Ω ── MOSFET gate
GND común ESP32 y fuente 12V (estrella)
```

**Nunca** alimentar solenoides desde el pin del ESP32.

---

## 3. Firmware

Archivo: `kilometro_v11_esp32.ino`

### 3.1 Comandos Serial (115200 baud)

| Cmd | Acción |
|-----|--------|
| `s` | START un ciclo (si IDLE_ALTA) |
| `r` | RESET a FAULT clear / rearmar |
| `p` | Print estado |
| `z` | Zero contadores de ciclo |
| `h` | Help |

### 3.2 Línea de log CSV

```text
cycle,t_enganche_ms,t_bajada_ms,t_entrega_ms,t_subida_ms,t_total_ms,ok,fault_code
```

Copiar del Monitor Serie a un `.csv` para análisis.

### 3.3 Timeouts (configurables)

| Parámetro | Default | Significado |
|-----------|---------|-------------|
| `HOLD_MS` | 300 | tiempo con ambos pernos cerrados |
| `SETTLE_MS` | 150 | asentamiento |
| `T_BAJADA_MAX_MS` | 30000 | timeout hundimiento |
| `T_SUBIDA_MAX_MS` | 45000 | timeout subida por boyancia |
| `SOL_PULSE_MAX_MS` | 2000 | no dejar solenoide en ON eterno (duty) |

Si el solenoide es de **servicio intermitente**, el firmware usa pulsos y opcionalmente “hold” PWM bajo (v1.1 básica: solo ON corto + pasador mecánico con muesca que se queda enclavado sin corriente — **preferible**).

### 3.4 Diseño mecánico preferido: enclavamiento biestable

Ideal para no calentar bobinas:

1. Solenoide **empuja** el pasador a posición ON o OFF.  
2. Un **muelle + muesca** mantiene la posición sin corriente.  
3. Firmware solo da **pulso** 100–400 ms.

Si el solenoide debe quedar energizado, respeta duty del fabricante y `SOL_PULSE_MAX_MS`.

---

## 4. Procedimiento de puesta en marcha

1. Montar mecánica v1 y validar 20 ciclos **manuales**.  
2. Instalar bowden + solenoides en seco; verificar carrera de cada pasador.  
3. Cablear MOSFET + flyback + GND común.  
4. Alimentar ESP32 por USB; abrir Serial 115200.  
5. `p` → estado IDLE_ALTA solo si FC_ALTA activo y módulo arriba.  
6. Sin lastre: probar enganche/entrega en aire (simular FC con la mano).  
7. Con agua y lastre: un ciclo con supervisión visual.  
8. 20 ciclos automáticos; exportar log.

---

## 5. Seguridad firmware

- No conmutar si el FC de cota no está activo (salvo modo `FORCE` de debug, desactivado por defecto).  
- Timeout → `FAULT`, LED rojo, no seguir.  
- Botón RESET limpia FAULT solo si el operador confirma posición segura.  
- Watchdog opcional (`enableLoopWDT`).

---

## 6. Qué medir en v1.1

| Magnitud | Cómo |
|----------|------|
| t_enganche, t_bajada, t_entrega, t_subida | log ESP32 |
| Fiabilidad transferencias | % ok / 50 ciclos |
| Temperatura solenoides | tacto / IR tras 20 ciclos |
| Corriente 12 V | pinza amperimétrica |

Balance energético del lastre sigue siendo el de la doc v1 (`m g Δh`); el ESP32 **no inventa julios**, solo **cronometra** el shuttle.

---

## 7. Archivos de esta carpeta

| Archivo | Descripción |
|---------|-------------|
| `README_v1_1.md` | Este documento |
| `kilometro_v11_esp32.ino` | Firmware ESP32 (Arduino) |
| `wiring_v1_1.txt` | Resumen de conexiones |
| `BOM_v1_1_extra.csv` | Electrónica adicional |

---

## 8. Siguiente (v1.2, no incluido)

- SD card log  
- Telemetría WiFi (Home Assistant)  
- Lift automático del lastre BAJA→ALTA  
- PWM hold en solenoides  
- Interlock óptico de presencia de lastre  

---

*v1.1 — solenoides + ESP32 — `kilometro_sim/v1_1_ESP32/`*
