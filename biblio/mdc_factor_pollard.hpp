/**
 * ============================================================================
 * MÒDUL: mdc_factor_pollard.hpp
 * Factorització per Pollard-Rho (variant Brent) amb mulmod64 segur.
 * Depèn de: aritmetica_modular.hpp, primalitat.hpp
 *
 * QUAN USAR:
 *   Complementa factorizar_mdc (k-sweep MDC) per als casos DESEQUILIBRATS:
 *   N = p·q on p << q (factor petit molt diferent del gran).
 *   El k-sweep és fort quan p ≈ q; Pollard-Brent és fort quan p << q.
 *
 * QUAN NO USAR:
 *   N > 2^64 (necessitaria BigInt), o quan factorizar_mdc ja troba el factor.
 *
 * COMPLEXITAT: O(N^(1/4)) esperada, probabilista.
 * ============================================================================
 */
#pragma once
#include <cstdint>
#include <cmath>
#include <numeric>
#include "aritmetica_modular.hpp"
#include "primalitat.hpp"

namespace mdc {

/**
 * pollard_brent_factor(N, c):
 *   Un sol intent amb la constant c en la funció f(x) = x²+c mod N.
 *   Retorna un factor no trivial de N, o 0 si falla (cicle sense factor).
 *
 *   Algoritme: Brent (1980), millora de Pollard-Rho en nombre de crides.
 *   f(x) = (x·x + c) mod N  via mulmod64 segur (sense overflow).
 */
static inline uint64_t pollard_brent_factor(uint64_t N, uint64_t c) {
    if (N % 2 == 0) return 2;

    uint64_t x = 2, y = 2, d = 1;
    uint64_t m = 128;   // mida del bloc de producte
    uint64_t q = 1;
    uint64_t ys = 0, r = 1;

    auto f = [&](uint64_t v) -> uint64_t {
        return (mulmod64(v, v, N) + c) % N;
    };

    while (d == 1) {
        x = y;
        for (uint64_t i = 0; i < r; ++i) y = f(y);
        uint64_t k = 0;
        while (k < r && d == 1) {
            ys = y;
            for (uint64_t i = 0; i < m && i < r - k; ++i) {
                y = f(y);
                uint64_t diff = (x > y) ? x - y : y - x;
                q = mulmod64(q, diff, N);
            }
            d = std::gcd(q, N);
            k += m;
        }
        r <<= 1;
        if (r > (1ULL << 30)) return 0;  // evita bucle infinit per a N molt grans
    }

    if (d == N) {
        // Cicle detectat — retrocedir pas a pas (fallback de Floyd)
        d = 1;
        y = ys;
        while (d == 1) {
            ys = f(ys);
            uint64_t diff = (x > ys) ? x - ys : ys - x;
            d = std::gcd(diff, N);
        }
    }

    return (d == N) ? 0 : d;
}

/**
 * pollard_brent(N):
 *   Intenta factoritzar N provant múltiples constants c = 1, 2, 3, ...
 *   Retorna el factor menut trobat, o 0 si N és primer o falla.
 *
 *   Primer comprova primalitat (Miller-Rabin) per no perdre temps.
 *   Prova fins a MAX_INTENTS constants c.
 */
static inline uint64_t pollard_brent(uint64_t N) {
    if (N < 4)      return 0;
    if (N % 2 == 0) return 2;
    if (is_prime(N)) return 0;  // N és primer, no factoritzable

    // Primers petits ràpids abans de Pollard (evita bucles per factors trivials)
    static const uint64_t PETITS[] = {3,5,7,11,13,17,19,23,29,31,37};
    for (uint64_t p : PETITS) {
        if (N % p == 0) return p;
    }

    const int MAX_INTENTS = 20;
    for (int c = 1; c <= MAX_INTENTS; ++c) {
        uint64_t f = pollard_brent_factor(N, static_cast<uint64_t>(c));
        if (f != 0 && f != N) {
            return std::min(f, N / f);
        }
    }
    return 0;  // tots els intents fallits
}

/**
 * factorizar_complet(N):
 *   Estratègia combinada: k-sweep MDC primer (determinista, fort per p≈q),
 *   Pollard-Brent si el k-sweep no troba res (fort per p<<q).
 *
 *   Per usar, inclou també factoritzacio_mdc.hpp.
 */
#ifdef MDC_FACTORITZACIO_INCLOSA
static inline uint64_t factorizar_complet(uint64_t N) {
    uint64_t f = factorizar_mdc(N);
    if (f) return f;
    return pollard_brent(N);
}
#endif

} // namespace mdc
