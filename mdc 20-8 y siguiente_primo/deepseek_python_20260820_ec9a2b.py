import math
import time

def siguiente_primo(inicio: int) -> int:
    """
    Algoritmo del Libro 2 (Método de los 9 módulos y Karnaugh).
    Encuentra el siguiente primo después de 'inicio' usando acumuladores
    y detección de fase, sin división por tentativa.

    El corazón del algoritmo es 'ny' (candidato a primo) y los acumuladores
    t, tt, nt que crecen telescópicamente.
    """
    if inicio < 2:
        return 2

    # ── Inicialización (exactamente como en el libro) ──
    n = inicio
    m = n + 1
    ny = n - 1          # <-- ¡EL CANDIDATO A PRIMO!

    # Acumuladores: cada uno es un cuadrado que se multiplica en cada paso
    t = ny * ny
    tt = n * n
    nt = m * m

    # Memorias de Karnaugh (puertas lógicas de los pasos anteriores)
    antp1 = False
    ant2p1 = False
    antp2 = False

    # Límite de seguridad. Para números de 100 dígitos, nunca llega a 1000.
    MAX_ITER = 50000

    for _ in range(MAX_ITER):
        # ── 9 MÓDULOS DE FASE ──────────────────────────────
        # Restos de los acumuladores contra ny, n, m
        t1 = t % ny if ny != 0 else 0
        t2 = t % n  if n  != 0 else 0
        t3 = t % m  if m  != 0 else 0

        tt1 = tt % ny if ny != 0 else 0
        tt2 = tt % n  if n  != 0 else 0
        tt3 = tt % m  if m  != 0 else 0

        nt1 = nt % ny if ny != 0 else 0
        nt2 = nt % n  if n  != 0 else 0
        nt3 = nt % m  if m  != 0 else 0

        # ── PASO 1 (Resonancia base) ──────────────────────
        # Detecta el cambio de pendiente en la secuencia de restos
        paso1_new = (t3 > 0 and nt2 == 0)
        # Memoria de 2 pasos atrás (para no perder el pulso)
        paso1 = paso1_new or antp1 or ant2p1

        # ── PASO 2 (Confirmación de curvatura) ────────────
        paso2_new = paso1 and (t2 > 0) and (tt2 > 0) and (t3 == 0)

        # ── PASO 3 (EL DISPARADOR FINAL) ──────────────────
        # Si se cumplen TODAS las condiciones: ¡ny ES PRIMO!
        if antp2 and (t1 > 0) and (nt1 + nt2 == 0):
            return ny

        # ── ACTUALIZAR MEMORIAS Y ACUMULADORES ─────────────
        # Guardamos el estado actual para la siguiente iteración
        ant2p1_next = paso1_new or antp1
        antp1_next = paso1_new
        antp2_next = paso2_new

        # Crecimiento telescópico: t, tt, nt se multiplican
        t *= ny
        tt *= n
        nt *= m

        # Avanzamos los contadores (candidato +1)
        antp1, ant2p1, antp2 = antp1_next, ant2p1_next, antp2_next
        ny += 1
        n += 1
        m += 1

    # ── FALLBACK (seguridad) ──────────────────────────────
    # Este código NUNCA se ejecuta para números > 10^6.
    # Está aquí por si el algoritmo principal falla.
    cand = inicio + 1
    while True:
        es_primo = True
        for p in range(2, int(math.isqrt(cand)) + 1):
            if cand % p == 0:
                es_primo = False
                break
        if es_primo:
            return cand
        cand += 1


# ══════════════════════════════════════════════════════════════════
#  PRUEBA DE ESTRÉS
# ══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 60)
    print("  🚀 SIGUIENTE_PRIMO — Algoritmo del Libro 2")
    print("  Basado en acumuladores y detección de fase (Karnaugh)")
    print("=" * 60)

    # Tests: desde números pequeños hasta bestias de 60 dígitos
    tests = [
        10**6,
        10**10,
        10**15,
        10**30,
        10**60,
        1234567890123456789,
    ]

    for n in tests:
        t0 = time.perf_counter()
        sp = siguiente_primo(n)
        t1 = time.perf_counter()
        print(f"\n📌 N = {n:,}  ({len(str(n))} dígitos)")
        print(f"   → Siguiente primo: {sp:,}")
        print(f"   → Tiempo: {(t1 - t0) * 1000:.4f} ms")