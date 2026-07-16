/**
 * ============================================================================
 * MÒDUL: mdc_crt.hpp
 * Teorema Xinès del Reste (CRT) i reconstrucció multimodular.
 * Depèn de: aritmetica_modular.hpp
 *
 * QUAN USAR:
 *   - Reconstruir un enter gran a partir dels seus residus mod m1, m2, ...
 *   - Combinar solucions diofàntiques per a mòduls coprimers.
 *   - Accelerar operacions aritmètiques grans reduint-les a mòduls petits.
 * ============================================================================
 */
#pragma once
#include <cstdint>
#include <vector>
#include <numeric>
#include "aritmetica_modular.hpp"

namespace mdc {

/**
 * crt_dos(r1,m1, r2,m2):
 *   Donats r1 ≡ x (mod m1) i r2 ≡ x (mod m2), amb gcd(m1,m2)=1,
 *   retorna x ≡ r (mod m1·m2) per l'algoritme de Garner.
 *   Retorna {residu, modul_combinat}, o {0,0} si els mòduls no son coprimers.
 */
struct ResultatCRT { uint64_t residu; uint64_t modul; };

static inline ResultatCRT crt_dos(uint64_t r1, uint64_t m1,
                                   uint64_t r2, uint64_t m2) {
    if (std::gcd(m1, m2) != 1) return {0, 0};
    int64_t inv = inverso_modular(static_cast<int64_t>(m1),
                                   static_cast<int64_t>(m2));
    if (inv < 0) return {0, 0};

    // x = r1 + m1 · ((r2 - r1) · inv_m1 mod m2)
    int64_t diff = (static_cast<int64_t>(r2) - static_cast<int64_t>(r1)) % static_cast<int64_t>(m2);
    if (diff < 0) diff += static_cast<int64_t>(m2);
    uint64_t t = mulmod64(static_cast<uint64_t>(diff),
                          static_cast<uint64_t>(inv), m2);
    uint64_t M = m1 * m2;  // assumim que no desborda uint64_t
    uint64_t x = (r1 + m1 * t) % M;
    return {x, M};
}

/**
 * crt_llista(residus, moduls):
 *   Generalitzat per a n mòduls coprimers dos a dos.
 *   Aplica crt_dos iterativament.
 *   Retorna {x, M} on M = producte de tots els mòduls.
 *   Si algun parell no és coprimer, retorna {0,0}.
 */
static inline ResultatCRT crt_llista(const std::vector<uint64_t>& residus,
                                      const std::vector<uint64_t>& moduls) {
    if (residus.size() != moduls.size() || residus.empty())
        return {0, 0};

    ResultatCRT acc = {residus[0], moduls[0]};
    for (size_t i = 1; i < residus.size(); ++i) {
        acc = crt_dos(acc.residu, acc.modul, residus[i], moduls[i]);
        if (acc.modul == 0) return {0, 0};
    }
    return acc;
}

/**
 * verificar_crt(x, residus, moduls):
 *   Comprova que x satisfà tots els residus donats.
 *   Útil per a tests: verifiquem la solució CRT sense recalcular.
 */
static inline bool verificar_crt(uint64_t x,
                                  const std::vector<uint64_t>& residus,
                                  const std::vector<uint64_t>& moduls) {
    if (residus.size() != moduls.size()) return false;
    for (size_t i = 0; i < residus.size(); ++i) {
        if (moduls[i] == 0) return false;
        if (x % moduls[i] != residus[i]) return false;
    }
    return true;
}

} // namespace mdc
