/**
 * ============================================================================
 * ReglaMecanicaUniversal.hpp — PROJECTE 33x1, v4 (modular)
 * Autor: Víctor Manzanares Alberola
 * Refactor: Claude (Anthropic) — juny 2026
 *
 * Esta classe ara és un AGREGADOR FI: tota la matemàtica viu en mòduls
 * independents sota mdc_lib/. La classe wrapper existeix només per
 * mantenir compatibilitat binària amb la DLL/VBA actual (mateixes
 * signatures de mètode que v3).
 *
 * Mòduls (cadascun compila i es pot provar sol):
 *   mdc_lib/aritmetica_modular.hpp  — mulmod64, powmod64, inverso_modular, raiz_entera
 *   mdc_lib/primalitat.hpp          — is_prime (Miller-Rabin determinista)
 *   mdc_lib/roda_p210.hpp           — criba per roda 2·3·5·7
 *   mdc_lib/factoritzacio_mdc.hpp   — k_sweep_mdc, factorizar_mdc
 *   mdc_lib/simetries_mdc.hpp       — es_punto_simetria, es_doble_modulo
 *   mdc_lib/regla_calculo.hpp       — ReglaCalculo (escala log + interpolació)
 *   mdc_lib/analisi_modular.hpp     — fase/curvatura/massa/rugositat/Goldbach/...
 * ============================================================================
 */
#pragma once

#include "mdc_lib/aritmetica_modular.hpp"
#include "mdc_lib/primalitat.hpp"
#include "mdc_lib/roda_p210.hpp"
#include "mdc_lib/factoritzacio_mdc.hpp"
#include "mdc_lib/simetries_mdc.hpp"
#include "mdc_lib/regla_calculo.hpp"
#include "mdc_lib/analisi_modular.hpp"

#include <cstdint>

// Reexporta l'estructura de solució diofàntica i l'espectre (compatibilitat).
using SolucionDiofantica = mdc::SolucionDiofantica;

struct EstadoConvergenciaAvanzada {
    double   d1, d2, d3, j_actual;
    uint32_t iteraciones;
};

struct EspectroResonancia {
    uint32_t cantidad_frecuencias;
    double   frecuencias[4];
};

/**
 * Classe wrapper: cada mètode delega directament a una funció del mòdul
 * corresponent. No conté lògica pròpia ni estat mutable (excepte la regla
 * de càlcul, que té la seua taula precomputada).
 */
class ReglaMecanicaUniversal {
private:
    mdc::ReglaCalculo regla_calculo_;

public:
    ReglaMecanicaUniversal(uint32_t precision, uint32_t tamano_escala)
        : regla_calculo_(precision, tamano_escala) {}

    // ── aritmetica_modular.hpp ──────────────────────────────────────────
    int64_t calcular_inverso_modular_cinematico(int64_t A, int64_t M) const {
        return mdc::inverso_modular(A, M);
    }
    uint64_t calcular_multiplicacion_modular_bits(uint64_t A, uint64_t B, uint64_t M) const {
        return (M == 0) ? 0 : mdc::mulmod64(A % M, B % M, M);
    }
    uint64_t calcular_raiz_cinematica_bits(uint64_t A) const {
        return mdc::raiz_entera(A);
    }

    // ── primalitat.hpp ───────────────────────────────────────────────────
    bool evaluar_ny_karnaugh(uint64_t ny) const {
        return mdc::is_prime(ny);
    }

    // ── factoritzacio_mdc.hpp ────────────────────────────────────────────
    uint64_t factorizar_mdc(uint64_t N) const {
        return mdc::factorizar_mdc(N);
    }

    // ── simetries_mdc.hpp ────────────────────────────────────────────────
    bool es_punto_simetria(uint64_t N, uint64_t p) const {
        return mdc::es_punto_simetria(N, p);
    }
    bool es_doble_modulo(uint64_t N, uint64_t D) const {
        return mdc::es_doble_modulo(N, D);
    }

    // ── regla_calculo.hpp ────────────────────────────────────────────────
    double multiplicar(uint32_t pos_a, uint32_t pos_b) const {
        return regla_calculo_.multiplicar(pos_a, pos_b);
    }
    double posicion_de(double x) const {
        return regla_calculo_.posicion_de(x);
    }
    double valor_en_posicion(double pos) const {
        return regla_calculo_.valor_en_posicion(pos);
    }
    double multiplicar_valores(double a, double b) const {
        return regla_calculo_.multiplicar_valores(a, b);
    }

    // ── analisi_modular.hpp ──────────────────────────────────────────────
    uint64_t buscar_siguiente_primo(uint64_t desde_ny) const {
        return mdc::buscar_siguiente_primo(desde_ny);
    }
    double calcular_fase_modular(uint64_t ny, uint64_t M) const {
        return mdc::calcular_fase_modular(ny, M);
    }
    uint64_t analizar_patron_bits(uint64_t ny, uint32_t ventana) const {
        return mdc::analizar_patron_bits(ny, ventana);
    }
    double calcular_curvatura_modular(uint64_t ny, uint64_t M) const {
        return mdc::calcular_curvatura_modular(ny, M);
    }
    double calcular_masa_modular(uint64_t ny, uint64_t M) const {
        return mdc::calcular_masa_modular(ny, M);
    }
    double analizar_rugosidad_multiescala(uint64_t ny, uint64_t M) const {
        return mdc::analizar_rugosidad_multiescala(ny, M);
    }
    bool detectar_inflexion_modular(uint64_t ny, uint64_t M) const {
        return mdc::detectar_inflexion_modular(ny, M);
    }
    bool escanear_goldbach_resonante(uint64_t N, uint64_t& p1, uint64_t& p2,
                                      uint64_t limite_barrido) const {
        return mdc::escanear_goldbach_resonante(N, p1, p2, limite_barrido);
    }
    bool detectar_primos_gemelos_fase(uint64_t ny, uint64_t& p1, uint64_t& p2) const {
        return mdc::detectar_primos_gemelos_fase(ny, p1, p2);
    }
    bool discriminant_exacte(uint64_t S, uint64_t N, uint64_t& p_out, uint64_t& q_out) const {
        return mdc::discriminant_exacte(S, N, p_out, q_out);
    }
    SolucionDiofantica resolver_diofantica_cinematico(
            int64_t A, int64_t B, int64_t C, uint64_t limite_barrido) const {
        return mdc::resolver_diofantica_cinematico(A, B, C, limite_barrido);
    }
};
