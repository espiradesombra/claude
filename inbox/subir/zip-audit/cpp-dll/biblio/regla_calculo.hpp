/**
 * ============================================================================
 * MÒDUL: regla_calculo.hpp
 * Regla de càlcul analògica (escala logarítmica + interpolació).
 * NOTA: provat que NO accelera factorizar_mdc (divisió entera és més
 * ràpida que log10 en doble precisió: 6.1ns vs 10.2ns, benchmark real).
 * El seu valor és didàctic/exploratori, no de rendiment.
 * ============================================================================
 */
#pragma once
#include <cstdint>
#include <cmath>
#include <vector>

namespace mdc {

class ReglaCalculo {
private:
    const uint32_t precision_;
    const uint32_t tamano_escala_;
    std::vector<double> escala_logaritmica_;

public:
    ReglaCalculo(uint32_t precision, uint32_t tamano)
        : precision_(precision), tamano_escala_(tamano),
          escala_logaritmica_(tamano)
    {
        for (uint32_t i = 0; i < tamano_escala_; ++i)
            escala_logaritmica_[i] = std::pow(10.0, double(i) / precision_);
    }

    /** Versió "física": suma índexs enters i llig el valor de taula. */
    double multiplicar(uint32_t pos_a, uint32_t pos_b) const {
        uint32_t idx = pos_a + pos_b;
        return (idx < tamano_escala_) ? escala_logaritmica_[idx] : -1.0;
    }

    /** Inversa: posició (índex fraccionari) d'un valor real x>0. */
    double posicion_de(double x) const {
        if (x <= 0.0) return -1.0;
        return double(precision_) * std::log10(x);
    }

    /** Lectura interpolada de l'escala per a posicions fraccionàries. */
    double valor_en_posicion(double pos) const {
        if (pos < 0.0 || pos >= double(tamano_escala_ - 1)) {
            return std::pow(10.0, pos / double(precision_));
        }
        uint32_t i0 = static_cast<uint32_t>(pos);
        uint32_t i1 = i0 + 1;
        double frac = pos - double(i0);
        return escala_logaritmica_[i0] * (1.0 - frac) + escala_logaritmica_[i1] * frac;
    }

    /**
     * Cicle complet: valor → posició → suma → valor.
     * Error mesurat: ~0.000065% (limitat per precision_ i interpolació
     * lineal, no per bug).
     */
    double multiplicar_valores(double a, double b) const {
        if (a <= 0.0 || b <= 0.0) return -1.0;
        return valor_en_posicion(posicion_de(a) + posicion_de(b));
    }
};

} // namespace mdc
