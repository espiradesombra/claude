/**
 * ============================================================================
 * MÒDUL: aritmetica_modular.hpp
 * Primitives d'aritmètica modular sense desbordament.
 * Sense estat — totes les funcions són pures (static inline).
 * ============================================================================
 */
#pragma once
#include <cstdint>

namespace mdc {

/**
 * mulmod64(a, b, m) = (a · b) mod m, sense desbordament.
 * Usa __uint128_t (GCC/Clang). Verificat amb operands prop del màxim
 * de uint64_t (vegeu banc de proves, TEST aritmètica).
 */
static inline uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<uint64_t>(
        (static_cast<__uint128_t>(a) * b) % m
    );
#else
    uint64_t res = 0;
    a %= m;
    while (b > 0) {
        if (b & 1) res = (res + a) % m;
        a = (a * 2 >= a) ? (a * 2 % m) : ((a - (m - a)) % m);
        b >>= 1;
    }
    return res;
#endif
}

/**
 * powmod64(base, exp, m) = base^exp mod m, exponentació ràpida O(log exp).
 */
static inline uint64_t powmod64(uint64_t base, uint64_t exp, uint64_t m) {
    if (m == 1) return 0;
    uint64_t res = 1;
    base %= m;
    while (exp > 0) {
        if (exp & 1) res = mulmod64(res, base, m);
        base = mulmod64(base, base, m);
        exp >>= 1;
    }
    return res;
}

/**
 * inverso_modular(A, M): Euclides estès, O(log M).
 * Retorna x tal que A·x ≡ 1 (mod M), o -1 si no existeix (gcd(A,M)≠1).
 * [FIX-5 original]: substitueix el bucle O(M) ingenu.
 */
static inline int64_t inverso_modular(int64_t A, int64_t M) {
    if (M <= 1) return -1;
    int64_t old_r = A % M, r = M;
    int64_t old_s = 1,     s = 0;
    while (r != 0) {
        int64_t q = old_r / r;
        int64_t tmp = r;   r   = old_r - q * r;   old_r = tmp;
        tmp = s;   s   = old_s - q * s;   old_s = tmp;
    }
    if (old_r != 1) return -1;
    return (old_s % M + M) % M;
}

/**
 * raiz_entera(A): arrel quadrada entera per Newton, exacta per a uint64_t.
 */
static inline uint64_t raiz_entera(uint64_t A) {
    if (A < 2) return A;
    uint64_t x = A >> 1;
    uint64_t prev = 0;
    while (x != prev) {
        prev = x;
        x = (x + A / x) >> 1;
    }
    return x;
}

} // namespace mdc
