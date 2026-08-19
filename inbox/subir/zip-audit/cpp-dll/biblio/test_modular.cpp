/**
 * test_modular.cpp — Banc de proves de la versió modular (v4).
 * Mateix contingut de tests que v3, ara contra els mòduls separats.
 */
#include "ReglaMecanicaUniversal.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <string>

int main() {
    ReglaMecanicaUniversal regla(1000, 50000);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "================================================================\n";
    std::cout << "  BANC DE PROVES — REGLA MECÀNICA UNIVERSAL v3\n";
    std::cout << "================================================================\n\n";

    // ── TEST 1: Primalitat Miller-Rabin ──────────────────────────────────
    std::cout << "[TEST 1] Miller-Rabin determinista:\n";
    struct { uint64_t n; bool esperat; } casos_mr[] = {
        {2,true},{3,true},{5,true},{7,true},{11,true},{13,true},
        {4,false},{9,false},{15,false},{100003,true},
        {100003ULL*100019ULL, false},
        {999999999999999989ULL, true}  // primer gran conegut
    };
    int ok_mr = 0;
    for (auto& c : casos_mr) {
        bool res = regla.evaluar_ny_karnaugh(c.n);
        bool pass = (res == c.esperat);
        if (pass) ++ok_mr;
        std::cout << "  " << (pass?"✓":"✗") << "  n=" << c.n
                  << "  primer=" << (res?"SI":"NO")
                  << (pass?"":" ← ERROR") << "\n";
    }
    std::cout << "  Encerts: " << ok_mr << "/" << sizeof(casos_mr)/sizeof(*casos_mr) << "\n\n";

    // ── TEST 2: Cerca de primers ──────────────────────────────────────────
    std::cout << "[TEST 2] Proper primer (roda p210):\n";
    struct { uint64_t des; uint64_t esp; } casos_bp[] = {
        {1200, 1201}, {1000, 1009}, {10000, 10007},
        {999999999999999900ULL, 999999999999999967ULL}
    };
    for (auto& c : casos_bp) {
        auto t0 = std::chrono::high_resolution_clock::now();
        uint64_t res = regla.buscar_siguiente_primo(c.des);
        auto t1 = std::chrono::high_resolution_clock::now();
        long us = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
        bool pass = (res == c.esp);
        std::cout << "  " << (pass?"✓":"✗")
                  << "  des=" << c.des << "  → " << res
                  << "  (" << us << " µs)"
                  << (pass?"":" ← ERROR (esperat: "+std::to_string(c.esp)+")")
                  << "\n";
    }
    std::cout << "\n";

    // ── TEST 3: K-sweep MDC — factorització ───────────────────────────────
    std::cout << "[TEST 3] factorizar_mdc (k-sweep MDC):\n";
    struct { uint64_t N; uint64_t p_esp; } casos_mdc[] = {
        {3ULL*5,                 3},
        {101ULL*103,             101},
        {100003ULL*100019ULL,    100003},
        {1000003ULL*1000033ULL,  1000003},
        {9999991ULL*9999973ULL,  9999973},
        {4294967291ULL*4294967279ULL, 4294967279ULL}
    };
    int ok_mdc = 0;
    for (auto& c : casos_mdc) {
        auto t0 = std::chrono::high_resolution_clock::now();
        uint64_t f = regla.factorizar_mdc(c.N);
        auto t1 = std::chrono::high_resolution_clock::now();
        long us = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
        bool pass = (f == c.p_esp);
        if (pass) ++ok_mdc;
        std::cout << "  " << (pass?"✓":"✗")
                  << "  N=" << c.N << "  → p=" << f
                  << "  (" << us << " µs)"
                  << (pass?"":" ← ERROR (esperat: "+std::to_string(c.p_esp)+")")
                  << "\n";
    }
    std::cout << "  Encerts: " << ok_mdc << "/" << sizeof(casos_mdc)/sizeof(*casos_mdc) << "\n\n";

    // ── TEST 4: Goldbach ──────────────────────────────────────────────────
    std::cout << "[TEST 4] Goldbach resonant:\n";
    for (uint64_t N : {uint64_t(100), uint64_t(840), uint64_t(1000), uint64_t(10000)}) {
        uint64_t p1=0, p2=0;
        bool ok = regla.escanear_goldbach_resonante(N, p1, p2, N/2);
        if (ok)
            std::cout << "  " << N << " = " << p1 << " + " << p2
                      << "  " << (p1+p2==N && regla.evaluar_ny_karnaugh(p1)
                                             && regla.evaluar_ny_karnaugh(p2) ? "✓":"✗")
                      << "\n";
        else
            std::cout << "  " << N << " → no trobat\n";
    }
    std::cout << "\n";

    // ── TEST 5: Invers modular ────────────────────────────────────────────
    std::cout << "[TEST 5] Invers modular (Euclides estès):\n";
    struct { int64_t A, M, esp; } casos_inv[] = {
        {3, 7, 5}, {17, 3120, 2753}, {7, 26, 15}, {4, 9, 7}
    };
    for (auto& c : casos_inv) {
        int64_t res = regla.calcular_inverso_modular_cinematico(c.A, c.M);
        bool pass = (res == c.esp);
        std::cout << "  " << (pass?"✓":"✗")
                  << "  " << c.A << "^-1 mod " << c.M << " = " << res
                  << (pass?"":" ← ERROR (esperat: "+std::to_string(c.esp)+")")
                  << "\n";
    }
    std::cout << "\n";

    // ── TEST 6: Primers bessons ───────────────────────────────────────────
    std::cout << "[TEST 6] Primers bessons (>= 1000):\n";
    uint64_t tw1=0, tw2=0;
    regla.detectar_primos_gemelos_fase(1000, tw1, tw2);
    std::cout << "  Parella: " << tw1 << ", " << tw2
              << (tw2 == tw1+2 && regla.evaluar_ny_karnaugh(tw1)
                               && regla.evaluar_ny_karnaugh(tw2) ? " ✓":" ✗")
              << "\n\n";

    std::cout << "================================================================\n";

    // ── TEST 7: Regla de càlcul analògica (multiplicar_valores) ──────────
    std::cout << "\n[TEST 7] Regla de càlcul (multiplicar_valores):\n";
    struct { double a, b; } casos_regla[] = {
        {2.0, 3.0}, {7.0, 8.0}, {12.5, 4.0}, {100.0, 100.0}, {3.14159, 2.71828}
    };
    double error_max = 0.0;
    for (auto& c : casos_regla) {
        double res = regla.multiplicar_valores(c.a, c.b);
        double real = c.a * c.b;
        double err_pct = std::abs(res - real) / real * 100.0;
        error_max = std::max(error_max, err_pct);
        std::cout << "  " << c.a << " × " << c.b << " = " << res
                  << "  (real: " << real << ", error: " << err_pct << "%)\n";
    }
    std::cout << "  Error màxim observat: " << error_max << "%\n";

    // ── TEST 8: Punt de simetria N mod 2p = p ─────────────────────────────
    std::cout << "\n[TEST 8] es_punto_simetria (N mod 2p = p):\n";
    struct { uint64_t N, p; bool esperat; } casos_sim[] = {
        {10403, 101, true},   // 101*103, p factor real
        {10403, 103, true},   // l'altre factor també
        {10403, 7,   false},  // 7 NO divideix 10403
        {10403, 50,  false},  // no és factor
        {77,    7,   true},
        {77,    11,  true},
    };
    int ok_sim = 0;
    for (auto& c : casos_sim) {
        bool res = regla.es_punto_simetria(c.N, c.p);
        bool pass = (res == c.esperat);
        if (pass) ++ok_sim;
        std::cout << "  " << (pass?"✓":"✗")
                  << "  N=" << c.N << " p=" << c.p
                  << " → simetria=" << (res?"SI":"NO")
                  << (pass?"":" ← ERROR") << "\n";
    }
    std::cout << "  Encerts: " << ok_sim << "/" << sizeof(casos_sim)/sizeof(*casos_sim) << "\n";

    // ── TEST 9: Identitat de doble mòdul N mod(N//D-1)=D ──────────────────
    std::cout << "\n[TEST 9] es_doble_modulo (asimètric, només factor petit):\n";
    struct { uint64_t N, D; bool esperat; } casos_dm[] = {
        {10403, 101, true},   // p (factor petit) -> SI
        {10403, 103, false},  // q (factor gran) -> NO (asimetria esperada)
        {10403, 100, false},  // no-factor
        {77, 7, true},
        {77, 11, false},
        {15, 3, true},
    };
    int ok_dm = 0;
    for (auto& c : casos_dm) {
        bool res = regla.es_doble_modulo(c.N, c.D);
        bool pass = (res == c.esperat);
        if (pass) ++ok_dm;
        std::cout << "  " << (pass?"✓":"✗")
                  << "  N=" << c.N << " D=" << c.D
                  << " → " << (res?"SI":"NO")
                  << (pass?"":" ← ERROR") << "\n";
    }
    std::cout << "  Encerts: " << ok_dm << "/" << sizeof(casos_dm)/sizeof(*casos_dm) << "\n";

    return 0;
}
