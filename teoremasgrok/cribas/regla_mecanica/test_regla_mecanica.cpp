#include <iostream>
#include "regla_mecanica_pendiente.hpp"
#include "ReglaMecanicaUniversal.hpp"

int main() {
    std::cout << "=== Regla mecanica (2v+3)(2l+3) ===\n";

    uint64_t x = 100;
    std::cout << "x=" << x << " choque_mitad=" << rm::choque_mitad(x) << "\n";
    for (uint32_t p = 2; p <= 5; ++p) {
        x = rm::encoger_en_choque(x, p);
        std::cout << "  tras pendiente " << p << " -> x=" << x << "\n";
    }

    struct {
        uint64_t N;
        uint64_t esp;
    } casos[] = {
        {15, 3},
        {101ULL * 103, 101},
        {100003ULL * 100019ULL, 100003},
    };

    ReglaMecanicaUniversal regla(1000, 50000);
    for (auto &c : casos) {
        uint64_t f_rm = rm::factorizar_por_pendiente(c.N);
        uint64_t f_mdc = regla.factorizar_mdc(c.N);
        std::cout << "N=" << c.N << " regla=" << f_rm << " mdc=" << f_mdc
                  << (f_rm == c.esp ? " OK" : " FAIL") << "\n";
    }
    return 0;
}