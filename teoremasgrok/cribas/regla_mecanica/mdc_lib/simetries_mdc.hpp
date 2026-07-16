/**
 * ============================================================================
 * MÒDUL: simetries_mdc.hpp
 * Identitats geomètriques del sawtooth accordion (descobriments de Víctor).
 * Sense dependències.
 * ============================================================================
 */
#pragma once
#include <cstdint>

namespace mdc {

/**
 * es_punto_simetria(N, p): N mod 2p = p.
 * Demostració: N mod 2p = p ⟺ N = p·(2k+1) ⟺ p|N i N/p senar.
 * Per a semiprims senars, equival exactament a "p és factor".
 * Detecta AMBDÓS factors (p i q) simètricament.
 * Verificat sense falsos positius (escombrat exhaustiu D∈[3,200) per N=10403).
 */
static inline bool es_punto_simetria(uint64_t N, uint64_t p) {
    if (p == 0) return false;
    return (N % (2 * p)) == p;
}

/**
 * es_doble_modulo(N, D): N mod (N//D − 1) = D.
 * ASIMÈTRIC: només detecta el factor PETIT p, no el gran q.
 * (Comprovat: N=10403=101·103, D=101→true, D=103→false.)
 * Útil com a verificació creuada independent específica del factor petit,
 * complementària a es_punto_simetria (no redundant).
 */
static inline bool es_doble_modulo(uint64_t N, uint64_t D) {
    if (D == 0) return false;
    uint64_t q_aprox = N / D;
    if (q_aprox == 0) return false;
    uint64_t divisor = q_aprox - 1;
    if (divisor == 0) return false;
    return (N % divisor) == D;
}

} // namespace mdc
