/**
 * ============================================================================
 * PROYECTO 33x1: REGLA MECÁNICA UNIVERSAL
 * NÚCLEO MATEMÁTICO HÍBRIDO (CINEMÁTICO, ARMÓNICO Y BINARIO)
 * Autor: Víctor
 * Archivo: ReglaMecanicaUniversal_Consolidado.cpp
 * ============================================================================
 */

#include <iostream>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <iomanip>

// ============================================================================
// 1. ESTRUCTURAS DE DATOS Y DEFINICIÓN DE LA CLASE
// ============================================================================

struct SolucionDiofantica {
    bool encontrada;
    int64_t x;
    int64_t y;
};

struct EstadoConvergenciaAvanzada {
    double d1;        // j (Estado anterior inmediato)
    double d2;        // jant (Estado penúltimo)
    double d3;        // Antepenúltimo estado
    double j_actual;  // Estimador en curso
    uint32_t iteraciones;
};

struct EspectroResonancia {
    uint32_t cantidad_frecuencias;
    double frecuencias[4]; 
};

class ReglaMecanicaUniversal {
private:
    const uint32_t precision;
    const uint32_t tamano_escala;
    double* escala_logaritmica;

    // Pipeline Temporal de Karnaugh (Memorias Analógicas)
    bool antp1;
    bool antp2;
    bool ant2p1;

    // Acumulador Maestro del Núcleo ny (Casi-Factoriales)
    uint64_t t_k;

public:
    ReglaMecanicaUniversal(uint32_t prec, uint32_t tamano);
    ~ReglaMecanicaUniversal();

    // Núcleo Geométrico Base
    void precargar_escala();
    double multiplicar(uint32_t pos_a, uint32_t pos_b);
    uint32_t calcular_freno_euler(uint64_t numero_entrada);
    
    // Motor Síncrono de Primalidad y Operadores de Bits
    bool evaluar_ny_karnaugh(uint64_t ny);
    uint64_t buscar_siguiente_primo(uint64_t desde_ny);
    uint64_t aplicar_regla_desplazamiento_origen(uint64_t ny, int modo);

    // Nuevas Funciones Avanzadas e Interferencia Modular
    SolucionDiofantica resolver_diofantica_cinematico(int64_t A, int64_t B, int64_t C, uint64_t limite_barrido);
    double integrar_segmento_regla(uint32_t pos_inicio, uint32_t pos_fin, double R_Ohms, double C_Faradios);
    bool filtrar_armonico_zypyzape(double frecuencia_red, double frecuencia_turbina);
    double calcular_log_inflexion(double a, double b, EstadoConvergenciaAvanzada& estado);
    EspectroResonancia calcular_autovalores_resonancia(double sistema[2][2], uint64_t limite_barrido);
    
    // Aritmética Avanzada por Desplazamiento Limpio
    uint64_t calcular_multiplicacion_modular_bits(uint64_t A, uint64_t B, uint64_t M);
    uint64_t calcular_raiz_cinematica_bits(uint64_t A);
};

// ============================================================================
// 2. IMPLEMENTACIÓN DE MÉTODOS DE LA CLASE
// ============================================================================

ReglaMecanicaUniversal::ReglaMecanicaUniversal(uint32_t prec, uint32_t tamano) 
    : precision(prec), tamano_escala(tamano), antp1(false), antp2(false), ant2p1(false), t_k(2) {
    escala_logaritmica = new double[tamano_escala];
    precargar_escala();
}

ReglaMecanicaUniversal::~ReglaMecanicaUniversal() {
    delete[] escala_logaritmica;
}

void ReglaMecanicaUniversal::precargar_escala() {
    for (uint32_t i = 0; i < tamano_escala; ++i) {
        escala_logaritmica[i] = std::pow(10.0, static_cast<double>(i) / precision);
    }
}

double ReglaMecanicaUniversal::multiplicar(uint32_t pos_a, uint32_t pos_b) {
    uint32_t indice_resultado = pos_a + pos_b;
    if (indice_resultado < tamano_escala) {
        return escala_logaritmica[indice_resultado];
    }
    return -1.0; 
}

uint32_t ReglaMecanicaUniversal::calcular_freno_euler(uint64_t numero_entrada) {
    const double e_const = 2.718281828459045;
    double raiz = std::sqrt(static_cast<double>(numero_entrada));
    return static_cast<uint32_t>(std::floor(3.0 - e_const * raiz));
}

bool ReglaMecanicaUniversal::evaluar_ny_karnaugh(uint64_t ny) {
    this->t_k *= ny;
    uint64_t residuo_t1 = this->t_k % ny;
    bool impacto_actual = (residuo_t1 == 0);

    bool es_compuesto = (impacto_actual && antp1) || (!antp2 && ant2p1);

    ant2p1 = antp2;
    antp2 = antp1;
    antp1 = impacto_actual;

    return !es_compuesto; 
}

uint64_t ReglaMecanicaUniversal::aplicar_regla_desplazamiento_origen(uint64_t ny, int modo) {
    uint64_t el_doble = ny << 1; 
    switch(modo) {
        case 1: return el_doble - ny;       
        case 2: return el_doble + el_doble; 
        case 3: return el_doble + ny;       
        default: return ny;
    }
}

uint64_t ReglaMecanicaUniversal::buscar_siguiente_primo(uint64_t desde_ny) {
    uint64_t candidato = (desde_ny % 2 == 0) ? desde_ny + 1 : desde_ny + 2;
    while (true) {
        if (evaluar_ny_karnaugh(candidato)) {
            return candidato;
        }
        candidato = aplicar_regla_desplazamiento_origen(candidato, 1) + 2; 
    }
}

SolucionDiofantica ReglaMecanicaUniversal::resolver_diofantica_cinematico(int64_t A, int64_t B, int64_t C, uint64_t limite_barrido) {
    SolucionDiofantica sol = {false, 0, 0};
    for (uint64_t ny = 1; ny <= limite_barrido; ++ny) {
        int64_t residuo = (C - (A * ny)) % B;
        if (residuo == 0) {
            sol.encontrada = true;
            sol.x = ny;
            sol.y = (C - (A * ny)) / B;
            break; 
        }
    }
    return sol;
}

double ReglaMecanicaUniversal::integrar_segmento_regla(uint32_t pos_inicio, uint32_t pos_fin, double R_Ohms, double C_Faradios) {
    if (pos_inicio >= pos_fin || pos_fin >= this->tamano_escala) return 0.0;
    double acumulador_voltaje = 0.0;
    double dt = 0.001; 

    for (uint32_t i = pos_inicio; i < pos_fin; ++i) {
        double v_inst = this->escala_logaritmica[i];
        acumulador_voltaje += (v_inst / R_Ohms) * (dt / C_Faradios);
    }
    return acumulador_voltaje; 
}

bool ReglaMecanicaUniversal::filtrar_armonico_zypyzape(double frecuencia_red, double frecuencia_turbina) {
    const double TOLERANCIA_FASE = 0.02;
    double desfase = std::fmod(frecuencia_red, frecuencia_turbina);
    bool colision_actual = (desfase < TOLERANCIA_FASE || desfase > (frecuencia_turbina - TOLERANCIA_FASE));
    
    bool armonico_detectado = (colision_actual && this->antp1) || (!this->antp2 && this->ant2p1);
    
    this->ant2p1 = this->antp2;
    this->antp2 = this->antp1;
    this->antp1 = colision_actual;
    
    return armonico_detectado;
}

double ReglaMecanicaUniversal::calcular_log_inflexion(double a, double b, EstadoConvergenciaAvanzada& estado) {
    if (std::abs(estado.j_actual - estado.d1) < 1e-7 || estado.iteraciones > 60) {
        return estado.j_actual;
    }
    estado.iteraciones++;

    estado.d3 = estado.d2;
    estado.d2 = estado.d1;
    estado.d1 = estado.j_actual;

    double b_elevado = std::pow(b, estado.j_actual + 1.0);
    double j_salto = (a - (a / b_elevado)) / a;
    estado.j_actual += j_salto;

    double denominador = estado.d2 - estado.d1;
    if (std::abs(denominador) > 1e-9) {
        double relacion_inflexion = (estado.d3 - estado.d2) / denominador;
        
        if (relacion_inflexion > 0.0) {
            double distancia_tramo = estado.j_actual - estado.d1;
            estado.j_actual += (2.0 * distancia_tramo); 
        }
    }
    return calcular_log_inflexion(a, b, estado);
}

EspectroResonancia ReglaMecanicaUniversal::calcular_autovalores_resonancia(double sistema[2][2], uint64_t limite_barrido) {
    EspectroResonancia espectro = {0, {0.0, 0.0, 0.0, 0.0}};
    const double TOLERANCIA_RESONANCIA = 0.005;
    double traza = sistema[0][0] + sistema[1][1];
    double determinante = (sistema[0][0] * sistema[1][1]) - (sistema[0][1] * sistema[1][0]);

    for (uint64_t ny = 1; ny <= limite_barrido; ++ny) {
        double lambda = static_cast<double>(ny) / 100.0; 
        double residuo_fase = (lambda * lambda) - (traza * lambda) + determinante;

        if (std::abs(residuo_fase) < TOLERANCIA_RESONANCIA) {
            if (espectro.cantidad_frecuencias == 0 || std::abs(espectro.frecuencias[espectro.cantidad_frecuencias - 1] - lambda) > 0.2) {
                espectro.frecuencias[espectro.cantidad_frecuencias] = lambda;
                espectro.cantidad_frecuencias++;
                if (espectro.cantidad_frecuencias >= 2) break; 
            }
        }
    }
    return espectro;
}

uint64_t ReglaMecanicaUniversal::calcular_multiplicacion_modular_bits(uint64_t A, uint64_t B, uint64_t M) {
    uint64_t acumulador_resultado = 0;
    A = A % M;

    while (B > 0) {
        if (B & 1) { 
            acumulador_resultado = (acumulador_resultado + A) % M;
        }
        A = (A << 1) % M; 
        B = B >> 1;
    }
    return acumulador_resultado;
}

uint64_t ReglaMecanicaUniversal::calcular_raiz_cinematica_bits(uint64_t A) {
    if (A == 0 || A == 1) return A;
    uint64_t estimador = A >> 1; 
    uint64_t d1_ant = 0;

    while (estimador != d1_ant) {
        d1_ant = estimador;
        estimador = (estimador + (A / estimador)) >> 1; 
    }
    return estimador;
}

// ============================================================================
// 3. BANCO DE PRUEBAS PRINCIPAL (MAIN TESTBENCH)
// ============================================================================

int main() {
    ReglaMecanicaUniversal regla(1000, 50000);
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "========================================================\n";
    std::cout << "   BANCO DE PRUEBAS: REGLA MECÁNICA UNIVERSAL (C++)\n";
    std::cout << "========================================================\n\n";

    // TEST 1: Primalidad síncrona
    std::cout << "[TEST 1] Buscando el siguiente primo desde 10000...\n";
    auto i1 = std::chrono::high_resolution_clock::now();
    uint64_t primo = regla.buscar_siguiente_primo(10000);
    auto f1 = std::chrono::high_resolution_clock::now();
    std::cout << "  > Primo hallado: " << primo << " en " 
              << std::chrono::duration_cast<std::chrono::microseconds>(f1 - i1).count() << " us.\n\n";

    // TEST 2: MDC Ecuación Diofántica
    std::cout << "[TEST 2] Resolviendo Ecuacion Diofantica (7x + 13y = 500)...\n";
    auto i2 = std::chrono::high_resolution_clock::now();
    SolucionDiofantica sol = regla.resolver_diofantica_cinematico(7, 13, 500, 1000);
    auto f2 = std::chrono::high_resolution_clock::now();
    if (sol.encontrada) {
        std::cout << "  > Coincidencia MDC: x = " << sol.x << ", y = " << sol.y << " (" 
                  << std::chrono::duration_cast<std::chrono::microseconds>(f2 - i2).count() << " us).\n\n";
    }

    // TEST 3: Logaritmo por Inflexión
    std::cout << "[TEST 3] Calculando log_2(1024) por aproximacion de inflexion...\n";
    EstadoConvergenciaAvanzada estado_log = {0.0, 0.0, 0.0, 1.0, 0};
    auto i3 = std::chrono::high_resolution_clock::now();
    double res_log = regla.calcular_log_inflexion(1024.0, 2.0, estado_log);
    auto f3 = std::chrono::high_resolution_clock::now();
    std::cout << "  > Resultado j: " << res_log << " (Iteraciones: " << estado_log.iteraciones << ") en "
              << std::chrono::duration_cast<std::chrono::microseconds>(f3 - i3).count() << " us.\n\n";

    // TEST 4: Resonancia y Autovalores
    std::cout << "[TEST 4] Escaneando frecuencias de resonancia (Matriz Estabilidad)...\n";
    double matriz[2][2] = {{3.0, 2.0}, {1.0, 4.0}};
    auto i4 = std::chrono::high_resolution_clock::now();
    EspectroResonancia esp = regla.calcular_autovalores_resonancia(matriz, 1000);
    auto f4 = std::chrono::high_resolution_clock::now();
    std::cout << "  > Frecuencias de acoplamiento: " << esp.cantidad_frecuencias << " (" 
              << std::chrono::duration_cast<std::chrono::microseconds>(f4 - i4).count() << " us).\n";
    for(uint32_t i=0; i<esp.cantidad_frecuencias; ++i) {
        std::cout << "    * Autovalor #" << (i+1) << ": " << esp.frecuencias[i] << " Hz\n";
    }
    std::cout << "\n";

    // TEST 5: Multiplicación Modular Rápida (Bits)
    std::cout << "[TEST 5] Multiplicacion Modular Rapida (Evitando Desborde)...\n";
    auto i5 = std::chrono::high_resolution_clock::now();
    uint64_t res_mod = regla.calcular_multiplicacion_modular_bits(7123456789ULL, 9876543210ULL, 1000000007ULL);
    auto f5 = std::chrono::high_resolution_clock::now();
    std::cout << "  > Resultado mod: " << res_mod << " en "
              << std::chrono::duration_cast<std::chrono::microseconds>(f5 - i5).count() << " us.\n\n";

    // TEST 6: Raíz Cinemática por Bits
    std::cout << "[TEST 6] Calculando Raiz Cuadrada Entera por bit-shifts...\n";
    auto i6 = std::chrono::high_resolution_clock::now();
    uint64_t res_raiz = regla.calcular_raiz_cinematica_bits(123456789012345ULL);
    auto f6 = std::chrono::high_resolution_clock::now();
    std::cout << "  > Resultado raiz entero: " << res_raiz << " en "
              << std::chrono::duration_cast<std::chrono::microseconds>(f6 - i6).count() << " us.\n\n";

    std::cout << "========================================================\n";
    return 0;
}
