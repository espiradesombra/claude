/**
 * ============================================================================
 * MÒDUL: primalitat.hpp
 * Test de primalitat determinista (Miller-Rabin) i cerca de primers.
 * Depèn de: aritmetica_modular.hpp
 * ============================================================================
 */
#pragma once
#include <cstdint>
#include "aritmetica_modular.hpp"

namespace mdc {

/**
 * miller_rabin_witness(n, a): true si 'a' és testimoni de composicitat.
 */
static inline bool miller_rabin_witness(uint64_t n, uint64_t a) {
    if (n % a == 0) return n != a;

    uint64_t d = n - 1;
    int r = 0;
    while ((d & 1) == 0) { d >>= 1; ++r; }

    uint64_t x = powmod64(a, d, n);
    if (x == 1 || x == n - 1) return false;

    for (int i = 0; i < r - 1; ++i) {
        x = mulmod64(x, x, n);
        if (x == n - 1) return false;
    }
    return true;
}

/**
 * is_prime(n): determinista, correcte per a tot n < 3.317·10^24
 * (cobreix tot uint64_t). Bases: {2,3,5,7,11,13,17,19,23,29,31,37}.
 * [FIX-1 original]: substitueix el pipeline Karnaugh incorrecte de v2.
 */
static inline bool is_prime(uint64_t n) {
    if (n < 2)  return false;
    if (n == 2 || n == 3 || n == 5 || n == 7 ||
        n == 11 || n == 13 || n == 17 || n == 19 ||
        n == 23 || n == 29 || n == 31 || n == 37) return true;
    if ((n & 1) == 0 || n % 3 == 0 || n % 5 == 0) return false;

    static const uint64_t bases[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (uint64_t a : bases) {
        if (miller_rabin_witness(n, a)) return false;
    }
    return true;
}

} // namespace mdc
