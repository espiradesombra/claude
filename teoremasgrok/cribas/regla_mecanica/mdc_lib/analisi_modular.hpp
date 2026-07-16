/**
 * ============================================================================
 * MÒDUL: analisi_modular.hpp
 * Funcions d'anàlisi modular: fase, curvatura, massa, rugositat, bessons,
 * Goldbach, discriminant exacte, diofàntica. Proves trivials o cap encara
 * (vegeu nota a cada funció) — usar amb cautela en producció.
 * Depèn de: primalitat.hpp, roda_p210.hpp
 * ============================================================================
 */
#pragma once
#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>
#include "primalitat.hpp"
#include "roda_p210.hpp"

namespace mdc {

/** fase(ny, M) = 2π·(ny mod M)/M. Blindada contra M=0. */
static inline double calcular_fase_modular(uint64_t ny, uint64_t M) {
    if (M == 0) return 0.0;
    constexpr double PI2 = 6.283185307179586;
    uint64_t residuo = ny % M;
    return PI2 * (double(residuo) / double(M));
}

/** Nombre de transicions de bit en una finestra de `ventana` bits de ny. */
static inline uint64_t analizar_patron_bits(uint64_t ny, uint32_t ventana) {
    if (ventana == 0 || ventana > 64) return 0;
    uint64_t transicions = 0;
    uint64_t mascara = (ventana < 64) ? ((1ULL << ventana) - 1) : ~0ULL;
    uint64_t bits = ny & mascara;
    bool ult = bits & 1;
    bits >>= 1;
    for (uint32_t i = 1; i < ventana; ++i) {
        bool act = bits & 1;
        if (act != ult) { ++transicions; ult = act; }
        bits >>= 1;
    }
    return transicions;
}

/** Curvatura discreta de la fase modular al voltant de ny. PROVA TRIVIAL. */
static inline double calcular_curvatura_modular(uint64_t ny, uint64_t M) {
    if (ny < 2 || M == 0) return 0.0;
    double f0 = calcular_fase_modular(ny - 1, M);
    double f1 = calcular_fase_modular(ny,     M);
    double f2 = calcular_fase_modular(ny + 1, M);
    double v1 = f1 - f0, v2 = f2 - f1;
    double acc = v2 - v1;
    double denom = std::pow(1.0 + v1 * v1, 1.5);
    return (std::abs(denom) < 1e-12) ? 0.0 : std::abs(acc) / denom;
}

/** massa(ny,M) = log10(M/(ny mod M)). PROVA TRIVIAL. */
static inline double calcular_masa_modular(uint64_t ny, uint64_t M) {
    if (M == 0) return 0.0;
    uint64_t r = ny % M;
    if (r == 0) return std::numeric_limits<double>::infinity();
    return std::log10(double(M) / double(r));
}

/** Combina patró de bits + fase + massa en una mètrica de rugositat. SENSE PROVA. */
static inline double analizar_rugosidad_multiescala(uint64_t ny, uint64_t M) {
    double r1 = double(analizar_patron_bits(ny, 8)) / 8.0;
    double fase = calcular_fase_modular(ny, M);
    double r2 = std::abs(std::sin(fase));
    double masa = calcular_masa_modular(ny, M);
    double r3 = std::isinf(masa) ? 1.0 : masa / (masa + 1.0);
    return std::sqrt(r1*r1 + r2*r2 + r3*r3);
}

/** Màxim local de curvatura. SENSE PROVA — llindar 0.01 arbitrari, no calibrat. */
static inline bool detectar_inflexion_modular(uint64_t ny, uint64_t M) {
    if (M == 0 || ny < 2) return false;
    double kappa_prev = calcular_curvatura_modular(ny - 1, M);
    double kappa_curr = calcular_curvatura_modular(ny,     M);
    double kappa_next = calcular_curvatura_modular(ny + 1, M);
    return (kappa_curr > kappa_prev) && (kappa_curr > kappa_next)
           && (kappa_curr > 0.01);
}

/** Cerca p1+p2=N (N parell). 4/4 casos provats. */
static inline bool escanear_goldbach_resonante(uint64_t N, uint64_t& p1, uint64_t& p2,
                                                uint64_t limite_barrido) {
    if (N % 2 != 0 || N <= 2) return false;
    uint64_t ny = 3;
    uint64_t lim = std::min(N / 2, limite_barrido);
    while (ny <= lim) {
        if (is_prime(ny)) {
            uint64_t comp = N - ny;
            if (is_prime(comp)) { p1 = ny; p2 = comp; return true; }
        }
        ny += (ny == 2) ? 1 : 2;
    }
    return false;
}

/** Cerca primers bessons p,p+2 amb p>=ny. Límit 10^6 passes. 1/1 provat. */
static inline bool detectar_primos_gemelos_fase(uint64_t ny, uint64_t& p1, uint64_t& p2) {
    if (ny < 3) ny = 3;
    uint64_t c = (ny % 2 == 0) ? ny + 1 : ny;
    for (uint64_t i = 0; i < 1'000'000; ++i) {
        if (is_prime(c) && is_prime(c + 2)) { p1 = c; p2 = c + 2; return true; }
        c = seguent_candidat_roda(c);
    }
    return false;
}

/** Buscar el proper primer des de desde_ny, amb roda p210. 4/4 casos provats. */
static inline uint64_t buscar_siguiente_primo(uint64_t desde_ny) {
    if (desde_ny <= 2) return 2;
    if (desde_ny <= 3) return 3;
    uint64_t c = (desde_ny % 2 == 0) ? desde_ny + 1 : desde_ny;
    for (uint64_t p : {uint64_t(3),uint64_t(5),uint64_t(7),uint64_t(11),uint64_t(13)}) {
        if (c == p) return p;
    }
    while (true) {
        if (is_prime(c)) return c;
        c = seguent_candidat_roda(c);
    }
}

struct SolucionDiofantica { bool encontrada; int64_t x; int64_t y; };

/** SENSE PROVA REAL — barrido lineal, O(limite_barrido). */
static inline SolucionDiofantica resolver_diofantica_cinematico(
        int64_t A, int64_t B, int64_t C, uint64_t limite_barrido) {
    SolucionDiofantica sol = {false, 0, 0};
    if (B == 0) return sol;
    for (uint64_t ny = 1; ny <= limite_barrido; ++ny) {
        int64_t residu = (C - A * static_cast<int64_t>(ny)) % B;
        if (residu == 0) {
            sol.encontrada = true;
            sol.x = static_cast<int64_t>(ny);
            sol.y = (C - A * sol.x) / B;
            break;
        }
    }
    return sol;
}

/** Mètode C (Libro 5): Δ(S)=S²+6S-(N-9)=k² → p,q. SENSE PROVA REAL. */
static inline bool discriminant_exacte(uint64_t S, uint64_t N,
                                        uint64_t& p_out, uint64_t& q_out) {
    if (N < 9 || S + 6 > N) return false;
    int64_t Ss = static_cast<int64_t>(S);
    int64_t delta = Ss*Ss + 6*Ss - static_cast<int64_t>(N - 9);
    if (delta < 0) return false;
    uint64_t k = static_cast<uint64_t>(std::sqrt(double(delta)));
    while (k > 0 && k * k > uint64_t(delta)) --k;
    while ((k+1)*(k+1) <= uint64_t(delta)) ++k;
    if (k * k != uint64_t(delta)) return false;

    int64_t Sp = Ss + static_cast<int64_t>(k);
    int64_t Sm = Ss - static_cast<int64_t>(k);
    if (Sp < 0 || Sm < 0) return false;
    if (Sp % 2 != 0 || Sm % 2 != 0) return false;

    uint64_t v = uint64_t(Sp) / 2;
    uint64_t b = uint64_t(Sm) / 2;
    uint64_t p = 2*v + 3, q = 2*b + 3;
    if (p > 1 && q > 1 && p != N && q != N && p * q == N) {
        p_out = std::min(p, q);
        q_out = std::max(p, q);
        return true;
    }
    return false;
}

} // namespace mdc
