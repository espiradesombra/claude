/**
 * dll_export.cpp — Export C para ReglaMecanicaUniversal.dll (Windows)
 * Compilar:
 *   g++ -shared -O2 -std=c++17 -DBUILDING_DLL dll_export.cpp -o ReglaMecanicaUniversal.dll -lm
 *   cl /LD /O2 /std:c++17 /DBUILDING_DLL dll_export.cpp /Fe:ReglaMecanicaUniversal.dll
 */
#define BUILDING_DLL 1

#include "ReglaMecanicaUniversal.hpp"
#include "regla_mecanica_pendiente.hpp"

#include <cstdint>

static ReglaMecanicaUniversal g_regla(1000, 50000);

extern "C" {

#ifdef _WIN32
#define RM_EXPORT __declspec(dllexport)
#else
#define RM_EXPORT __attribute__((visibility("default")))
#endif

RM_EXPORT uint64_t rm_factorizar_mdc(uint64_t N) {
    return g_regla.factorizar_mdc(N);
}

RM_EXPORT uint64_t rm_factorizar_regla_mecanica(uint64_t N) {
    return rm::factorizar_por_pendiente(N);
}

RM_EXPORT int rm_es_primo(uint64_t n) {
    return g_regla.evaluar_ny_karnaugh(n) ? 1 : 0;
}

RM_EXPORT uint64_t rm_siguiente_primo(uint64_t desde) {
    return g_regla.buscar_siguiente_primo(desde);
}

RM_EXPORT uint64_t rm_encoger(uint64_t x, uint32_t pendiente) {
    return rm::encoger_en_choque(x, pendiente);
}

RM_EXPORT uint64_t rm_choque_mitad(uint64_t x) {
    return rm::choque_mitad(x);
}

RM_EXPORT uint64_t rm_pendiente_desde_v(uint64_t N, uint64_t v) {
    return rm::pendiente_desde_v(N, v);
}

RM_EXPORT int rm_factor_par(uint64_t N, uint64_t *f1, uint64_t *f2) {
    rm::FactorPar p = rm::factorizar_par(N);
    if (!p.ok || !f1 || !f2) {
        return 0;
    }
    *f1 = p.f1;
    *f2 = p.f2;
    return 1;
}

} // extern "C"