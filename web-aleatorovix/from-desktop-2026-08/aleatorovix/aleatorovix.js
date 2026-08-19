/**
 * ALEATOROVIX v1.0 (navegador)
 * Port de nucleo/entropia.py + mascara_lila.py + ciclo del organismo.
 * Víctor Manzanares Alberola · VMA / 33×1
 *
 * SIN Math.random.
 * Flujo: entropía → máscara Lila → intérprete mutante → bit flor / u / int
 *
 * Uso:
 *   <script src="aleatorovix/aleatorovix.js"></script>
 *   Aleatorovix.bit()          // 0 | 1
 *   Aleatorovix.random()       // [0, 1)
 *   Aleatorovix.randomInt(n)   // [0, n)
 *   Aleatorovix.ciclo()        // { decision, medida, k, bit, u, m }
 */
(function (global) {
  'use strict';

  var ciclos = 0;
  // buffer “código de operar” que se borra tras cada ciclo (criba desmemoriada)
  var codigoOperar = 'Lila v1.0 - 10^{1/x} + gemelos + pila + mierda + olvido';

  function getNanos() {
    if (typeof performance !== 'undefined' && performance.now) {
      var t = performance.now();
      return Math.floor((t % 1000) * 1e6 + (t * 1e3) % 1e3) % 1000000000;
    }
    return ((Date.now() % 1000) * 1000000) % 1000000000;
  }

  function robarBitPila() {
    // en browser no hay stack pointer: reloj + identidad de un objeto local
    var local = { t: typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now() };
    var h = 0;
    var s = String(local.t) + String(ciclos);
    for (var i = 0; i < s.length; i++) {
      h = ((h << 5) - h + s.charCodeAt(i)) | 0;
    }
    try {
      h ^= String(local).length << 3;
    } catch (e) { /* ignore */ }
    return (h >> 6) & 1;
  }

  function robarBasuraMemoria() {
    var basura = 0;
    if (typeof crypto !== 'undefined' && crypto.getRandomValues) {
      // ruido de plataforma (CSPRNG del SO; no es Math.random)
      var a = new Uint32Array(1);
      crypto.getRandomValues(a);
      basura = a[0] >>> 0;
    } else {
      basura = (getNanos() ^ ((Date.now() * 2654435761) >>> 0)) >>> 0;
    }
    var punteroSucio = (basura ^ getNanos()) >>> 0;
    return punteroSucio & 0xFF;
  }

  /** jitter sintético: micro-trabajo en CPU (sustituto de sleep + ping) */
  function jitterLocal(inerciaA) {
    var n = 30 + (inerciaA % 70);
    var acc = 0;
    var t0 = typeof performance !== 'undefined' ? performance.now() : Date.now();
    for (var i = 0; i < n; i++) {
      acc = (acc + i * 17 + (inerciaA & 7)) | 0;
    }
    var t1 = typeof performance !== 'undefined' ? performance.now() : Date.now();
    var us = Math.max(1, Math.floor((t1 - t0) * 1000) || (acc & 0xFF));
    return us;
  }

  /** Señal de red si el navegador la expone (opcional). No bloquea si no hay red. */
  function senalRedOpcional() {
    try {
      var c = typeof navigator !== 'undefined' && (navigator.connection || navigator.mozConnection || navigator.webkitConnection);
      if (!c) return 0;
      var h = 0;
      if (typeof c.downlink === 'number') h ^= (c.downlink * 1000) | 0;
      if (typeof c.rtt === 'number') h ^= (c.rtt * 17) | 0;
      if (c.effectiveType) {
        var s = String(c.effectiveType);
        for (var i = 0; i < s.length; i++) h = ((h << 5) - h + s.charCodeAt(i)) | 0;
      }
      if (typeof navigator.onLine === 'boolean') h ^= navigator.onLine ? 0x5a5a : 0xa5a5;
      return h >>> 0;
    } catch (e) {
      return 0;
    }
  }

  function muestrear() {
    var nanos = getNanos();
    var inerciaA = nanos % 1000;
    var pingUs = jitterLocal(inerciaA);
    // red_x = jitter local + (si hay) huella de Network Information API / online
    var redX = ((pingUs % 1000) ^ (senalRedOpcional() % 997)) % 1000;
    return {
      nanos: nanos,
      inercia_a: inerciaA,
      red_x: redX || 1,
      bit_pila: robarBitPila(),
      basura: robarBasuraMemoria(),
      ping_us: pingUs
    };
  }

  /** (-10^(1/x) + 10^(1/(x+a)) - 1) * x  — organismo Lila */
  function mascaraLila(x, a) {
    if (x <= 0) x = 1;
    if (a < 0) a = 0;
    var t1 = -Math.pow(10, 1 / x);
    var t2 = Math.pow(10, 1 / (x + a));
    return (t1 + t2 - 1) * x;
  }

  function pulsoGauss(valor, pico) {
    return Math.abs(valor - pico) <= 1 ? 1 : 0;
  }

  function interpreteMutante(bitExterno, bitRobado, ruidoInercia) {
    var par = (bitExterno << 1) | (bitRobado & 1);
    var rotacion = ruidoInercia % 3;
    var base;
    if (par === 0b11 || par === 0b00) base = 0;
    else if (par === 0b01) base = 1;
    else base = 2; // 0b10
    return (base + rotacion) % 3;
  }

  function decidir(m) {
    var curva = mascaraLila(m.red_x, m.inercia_a);
    var medida = Math.abs(Math.floor(curva)) % 10;
    var b0 = pulsoGauss(medida, 0);
    var r0 = interpreteMutante(b0, m.bit_pila, m.basura);
    var b5 = pulsoGauss(medida, 5);
    var r5 = interpreteMutante(b5, m.bit_pila, m.basura ^ 1);
    var b9 = pulsoGauss(medida, 9);
    var r9 = interpreteMutante(b9, m.bit_pila, m.basura ^ (m.nanos & 0xFF));
    var decision = (r0 + r5 + r9) % 4;
    return {
      decision: decision,
      medida: medida,
      r0: r0,
      r5: r5,
      r9: r9,
      curva: curva,
      inercia_a: m.inercia_a,
      red_x: m.red_x
    };
  }

  function desmemoriar() {
    var len = codigoOperar.length || 8;
    codigoOperar = Array(len + 1).join('\0');
  }

  /**
   * Un ciclo completo del organismo (rtl = pinza desde arriba).
   * @returns {{ decision: number, medida: number, k: number, bit: number, u: number, m: object }}
   */
  function ciclo() {
    codigoOperar = 'Lila v1.0 - 10^{1/x} + gemelos + pila + mierda + olvido';
    var m = muestrear();
    var d = decidir(m);
    var kRaw = (
      m.nanos
      ^ (m.basura << 7)
      ^ (m.bit_pila << 15)
      ^ (m.inercia_a << 22)
      ^ ciclos
    ) % 10000;
    var k = 9999 - kRaw; // rtl
    var bit = (d.decision ^ m.bit_pila ^ (m.basura & 1) ^ (k & 1)) & 1;
    var u = (
      ((m.nanos ^ (m.basura * 0x9e3779b9) ^ (k * 0x85ebca6b) ^ (d.decision << 28)) >>> 0)
      / 4294967296
    );
    ciclos++;
    desmemoriar();
    return {
      decision: d.decision,
      medida: d.medida,
      k: k,
      bit: bit,
      u: u,
      m: m
    };
  }

  /** 0 o 1 — una flor */
  function bit() {
    return ciclo().bit;
  }

  /** drop-in de Math.random() pero de un ciclo Aleatorovix → [0, 1) */
  function random() {
    return ciclo().u;
  }

  /** entero en [0, max) */
  function randomInt(max) {
    if (max <= 0) return 0;
    if (max === 1) return 0;
    var c1 = ciclo();
    var c2 = ciclo();
    var acc = (
      (c1.k * 10007 + c2.k)
      ^ (c1.decision << 20)
      ^ (c2.medida << 12)
      ^ (c1.m.nanos & 0xFFFF)
    ) >>> 0;
    return acc % max;
  }

  function getCiclos() {
    return ciclos;
  }

  function reset() {
    ciclos = 0;
    codigoOperar = 'Lila v1.0 - 10^{1/x} + gemelos + pila + mierda + olvido';
  }

  var Aleatorovix = {
    version: '1.0-browser',
    ciclo: ciclo,
    bit: bit,
    random: random,
    randomInt: randomInt,
    getCiclos: getCiclos,
    reset: reset,
    mascaraLila: mascaraLila,
    decidir: decidir,
    muestrear: muestrear
  };

  // global (script tag) + CommonJS / AMD light
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = Aleatorovix;
  }
  global.Aleatorovix = Aleatorovix;
})(typeof globalThis !== 'undefined' ? globalThis : typeof window !== 'undefined' ? window : this);
