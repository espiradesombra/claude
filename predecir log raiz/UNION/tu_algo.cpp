#include "tu_algo.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

static inline double clamp(double x, double lo, double hi) {
    return std::min(std::max(x, lo), hi);
}

// Log_b(a): doble ajuste simultáneo con parámetros optimizados (FACTOR5 y FACTOR6 incluidos)
double tu_log_a_b(double a, double b) {
    if (!(a > 0.0) || !(b > 0.0) || b == 1.0) return std::numeric_limits<double>::quiet_NaN();

    // PARÁMETROS OPTIMIZADOS FINALES (de junto.docx)
    const double factor1 = 1.2;  // sobrepaso
    const double factor2 = 3.4;  // multiplicativo
    const double factor3 = -0.8; // exponencial
    const double factor4 = 1.9;  // exponencial
    const double factor5 = -0.8; // secundario (exponencial)
    const double factor6 = 1.9;  // secundario (exponencial)

    double j = 1.0;
    double d3 = 0.0, d2 = 0.0, d1 = 0.0; // Diferencias para el ajuste
    const double precision = 1e-15;
    const double limite_inferior = 1e-12;

    for (int iteracion = 0; iteracion < 15; ++iteracion) {
        double potencia = std::pow(b, j);
        // 1. Prevención de sobrepaso (overflow)
        if (a / potencia >= potencia) {
            j *= factor1;
            continue;
        }

        // Almacenar el valor anterior de j antes de la actualización principal
        double j_previo = j;

        // 2. CÁLCULO PRINCIPAL (Ecuación base)
        j = (a + (a / std::pow(b, j)) - 1.0) / a;

        // 3. ACTUALIZAR SECUENCIA DE DIFERENCIAS (d3, d2, d1)
        d3 = d2;
        d2 = d1;
        d1 = j - j_previo; // Diferencia actual

        // --- PRIMER AJUSTE (Multiplicativo) ---
        j *= (1.0 + factor2 * std::abs(d1));

        // --- SEGUNDO AJUSTE (Exponencial, condicional) ---
        if (std::abs(d3 - d1) > limite_inferior) {
            double factor_exp1 = factor3 + factor4 * (d3 - d2) / (d3 - d1);
            j = std::pow(std::abs(j), std::abs(factor_exp1));
        }

        // 4. SEGUNDA ACTUALIZACIÓN (Preparación para el tercer ajuste)
        j_previo = j;
        // Re-calcular j usando la fórmula base con el valor ajustado
        j = (a + (a / std::pow(b, j)) - 1.0) / a;

        // Actualizar las diferencias para el tercer ajuste
        d3 = d2;
        d2 = d1;
        d1 = j - j_previo;

        // --- TERCER AJUSTE (Exponencial, con factor5 y factor6 - Incondicional en el doc?) ---
        // NOTA: El documento muestra este ajuste sin condición 'if', pero añadimos una para evitar división por cero.
        if (std::abs(d3 - d1) > limite_inferior) {
            double factor_exp2 = factor5 + factor6 * (d3 - d2) / (d3 - d1);
            j = std::pow(std::abs(j), std::abs(factor_exp2));
        }

        // 5. COMPROBAR CONVERGENCIA
        if (std::abs(j - j_previo) < precision) {
            break;
        }
    }
    return j;
}


// Raíz n-ésima: x^n ≈ a, usando el método de doble ajuste simultáneo
double tu_root(double a, double n) {
    if (!(a >= 0.0) || !(n > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    if (a == 0.0) return 0.0;

    // PARÁMETROS OPTIMIZADOS FINALES (los mismos que para logaritmos)
    const double factor1 = 1.2;  // sobrepaso
    const double factor2 = 3.4;  // multiplicativo
    const double factor3 = -0.8; // exponencial
    const double factor4 = 1.9;  // exponencial
    const double factor5 = -0.8; // secundario (exponencial)
    const double factor6 = 1.9;  // secundario (exponencial)

    double x = std::max(std::pow(a, 1.0/n), 1e-300);
    double d3 = 0.0, d2 = 0.0, d1 = 0.0;
    const double precision = 1e-15;
    const double limite_inferior = 1e-12;

    for (int iteracion = 0; iteracion < 20; ++iteracion) {
        double potencia = std::pow(x, n);
        
        // 1. Prevención de sobrepaso (overflow)
        if (a / potencia >= potencia) {
            x *= factor1;
            continue;
        }

        // Almacenar el valor anterior
        double x_previo = x;

        // 2. CÁLCULO PRINCIPAL (Ecuación base adaptada para raíces)
        x = (a + (a / std::pow(x, n)) - 1.0) / a;

        // 3. ACTUALIZAR SECUENCIA DE DIFERENCIAS
        d3 = d2;
        d2 = d1;
        d1 = x - x_previo;

        // --- PRIMER AJUSTE (Multiplicativo) ---
        x *= (1.0 + factor2 * std::abs(d1));

        // --- SEGUNDO AJUSTE (Exponencial, condicional) ---
        if (std::abs(d3 - d1) > limite_inferior) {
            double factor_exp1 = factor3 + factor4 * (d3 - d2) / (d3 - d1);
            x = std::pow(std::abs(x), std::abs(factor_exp1));
        }

        // 4. SEGUNDA ACTUALIZACIÓN (Preparación para el tercer ajuste)
        x_previo = x;
        x = (a + (a / std::pow(x, n)) - 1.0) / a;

        // Actualizar diferencias para el tercer ajuste
        d3 = d2;
        d2 = d1;
        d1 = x - x_previo;

        // --- TERCER AJUSTE (Exponencial, con factor5 y factor6) ---
        if (std::abs(d3 - d1) > limite_inferior) {
            double factor_exp2 = factor5 + factor6 * (d3 - d2) / (d3 - d1);
            x = std::pow(std::abs(x), std::abs(factor_exp2));
        }

        // 5. COMPROBAR CONVERGENCIA
        if (std::abs(x - x_previo) < precision * std::max(1.0, std::abs(x_previo))) {
            break;
        }

        // Limitar valores extremos
        x = clamp(x, 1e-300, 1e+300);
    }
    return x;
}