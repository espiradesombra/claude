/**
 * ============================================================================
 * MÒDUL: roda_p210.hpp
 * Roda 2·3·5·7=210 per a iteració eficient de candidats coprimers.
 * Elimina el 77.1% de candidats sense cap divisió.
 * Sense dependències.
 * ============================================================================
 */
#pragma once
#include <cstdint>

namespace mdc {

/**
 * Salts entre residus consecutius vàlids mod 210 (coprimers amb 2,3,5,7).
 * 48 valors = φ(210). Verificat independentment amb generació per gcd.
 */
static const int SALTS_RODA_210[48] = {
    10,2,4,2,4,6,2,6,4,2,4,6,6,2,6,4,2,6,4,6,8,4,2,4,
    2,4,8,6,4,6,2,4,6,2,6,6,4,2,4,6,2,6,4,2,4,2,10,2
};

/**
 * Residus vàlids mod 210, en ordre. RESIDUS_210[1]==11 és l'inici típic
 * per a cribar des de D=11 (els primers 2,3,5,7 es tracten a banda).
 */
static const int RESIDUS_210[48] = {
    1,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,
    71,73,79,83,89,97,101,103,107,109,113,121,127,131,
    137,139,143,149,151,157,163,167,169,173,179,181,
    187,191,193,197,199,209
};

/**
 * seguent_candidat_roda(d): donat un senar d>=11, retorna el proper
 * candidat coprimer amb 210 (saltant compostos evidents).
 */
static inline uint64_t seguent_candidat_roda(uint64_t d) {
    int rem = static_cast<int>(d % 210);
    for (int i = 0; i < 48; ++i) {
        if (RESIDUS_210[i] == rem) return d + SALTS_RODA_210[i];
    }
    for (int delta = 2; delta < 210; delta += 2) {
        uint64_t nd = d + delta;
        int nr = static_cast<int>(nd % 210);
        for (int r : RESIDUS_210) if (r == nr) return nd;
    }
    return d + 2; // fallback, no hauria d'arribar ací
}

} // namespace mdc
