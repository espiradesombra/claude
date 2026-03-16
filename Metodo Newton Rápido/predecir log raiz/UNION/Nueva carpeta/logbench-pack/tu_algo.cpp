
#include "tu_algo.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

static inline double clamp(double x, double lo, double hi) {
    return std::min(std::max(x, lo), hi);
}

double tu_log_b_a(double a, double b) {
    if (!(a > 0.0) || !(b > 0.0) || b == 1.0) return std::numeric_limits<double>::quiet_NaN();
    if (a == 1.0) return 0.0; // Caso especial

    const double factor1 = 1.3;
    const double factor2 = 3.4;
    const double factor3 = -0.8;
    const double factor4 = 1.9;
    const double factor5 = -0.8;
    const double factor6 = 1.9;

    double j = 1.0;
    double d3 = 0.0, d2 = 0.0, d1 = 0.0;
    const double precision = 0.001;
    const double limite_inferior = 0.001;
    double potencia = std::pow(b, j);
    bool condicion = true;
    
    while (condicion == true) { // Aumenté iteraciones
        

        

        if (a / potencia >= potencia) {
            j++;
            j *= factor1;
            continue;
        }

        double j_previo = j;
        // Fórmula principal CORREGIDA
        j *= (a + (a / potencia) - 1.0) / a;

        d3 = d2;
        d2 = d1;
        d1 = j - j_previo;

        // Ajuste multiplicativo
        j *= (1.0 + factor2 * d1);

        // Ajuste exponencial 1
        
        double factor_exp1 = factor3 + factor4 * (d3 - d2) / (d3 - d1);
        j = std::pow((j), (factor_exp1));
        

        
        potencia = std::pow(b, j);
        j_previo = j;
        
        j *= (a + (a / potencia) - 1.0) / a;

        d3 = d2;
        d2 = d1;
        d1 = j - j_previo;

        // Ajuste exponencial 2
       
        double factor_exp2 = factor5 + factor6 * (d3 - d2) / (d3 - d1);
        j = std::pow((j), (factor_exp2));
              
        

        if (d1 < precision) { 
            condicion = false;
        }
        j = clamp(j, 1e-300, 1e300); // Evitar overflow/underflow
      
    }
    return j;
}

double tu_root(double a, double n) {
    if (!(a >= 0.0) || !(n > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    if (a == 0.0) return 0.0;
    if (a == 1.0) return 1.0;

    const double factor1 = 1.2;
    const double factor2 = 3.4;
    const double factor3 = -0.8;
    const double factor4 = 1.9;
    const double factor5 = -0.8;
    const double factor6 = 1.9;

    // Mejor inicialización
    double x = (a > 1.0) ? a : 1.0;
    double d3 = 0.0, d2 = 0.0, d1 = 0.0;
    const double precision = 0.001;
    const double limite_inferior = 0.001;
    bool condicion = true;
    while (condicion == true) {
        double potencia = std::pow(x, n);
        
        

        double x_previo = x;
        
        // FÓRMULA CORREGIDA - multiplicar por x
        x = x * (a + (a / potencia) - 1.0) / a;

        d3 = d2;
        d2 = d1;
        d1 = x - x_previo;

        // Ajustes
        x *= (1.0 + factor2 * std::abs(d1));

        if (std::abs(d3 - d1) > limite_inferior && d3 != d1) {
            double factor_exp1 = factor3 + factor4 * (d3 - d2) / (d3 - d1);
            x = std::pow((x), (factor_exp1));
        }

        x_previo = x;
        potencia = std::pow(x, n);
        if (std::isinf(potencia) || potencia == 0.0) break;
        
        x = x * (a + (a / potencia) - 1.0) / a;

        d3 = d2;
        d2 = d1;
        d1 = x - x_previo;

        if (std::abs(d3 - d1) > limite_inferior && d3 != d1) {
            double factor_exp2 = factor5 + factor6 * (d3 - d2) / (d3 - d1);
            x = std::pow((x), (factor_exp2));
        }

        if ( d1 < precision ) {
            condicion = false;
        }
        x = clamp(x,1e-300, 1e300); // Evitar overflow/underflow
        
        
    }
    return x;
}