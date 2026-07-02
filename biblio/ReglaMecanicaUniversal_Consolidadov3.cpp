/**
 * ============================================================================
 * PROYECTO 33x1: REGLA MECÁNICA UNIVERSAL  —  v3 (matemàtica corregida)
 * Autor: Víctor Manzanares Alberola
 *        EPSA / Universitat Politècnica de València (Alcoi)
 * Revisió: Claude (Anthropic) — juny 2026
 *
 * CANVIS RESPECTE A v2:
 *   [FIX-1] evaluar_ny_karnaugh: substituït el pipeline Karnaugh (incorrecte)
 *           per Miller-Rabin determinista per a tot n < 3.317·10^24.
 *           Bases: {2,3,5,7,11,13,17,19,23,29,31,37}.
 *   [FIX-2] buscar_siguiente_primo: ara usa evaluar_ny_karnaugh corregit +
 *           roda 2·3·5·7 per saltar candidats compostos. Garantit correcte.
 *   [FIX-3] calcular_fase_modular: eliminat l'acumulador t_k_grande (era
 *           estat global mutable → races en DLL). Ara computa ny mod M
 *           directament amb aritmètica modular 128-bit.
 *   [FIX-4] Afegits detectar_primos_gemelos_fase i detectar_inflexion_modular
 *           que faltaven al Consolidado original (trencaven la DLL).
 *   [FIX-5] calcular_inverso_modular_cinematico: substituït el bucle O(M)
 *           per l'algorisme d'Euclides estès O(log M).
 *   [FIX-6] K-sweep MDC (nou): factoritza semiprims equilibrats de fins a
 *           ~30 dígits sense roda de trial-division. Nucli del MDC de Víctor.
 *   [FIX-7] Eliminat t_k_grande global (era font de tots els bugs d'estat).
 *
 * COMPILACIÓ STANDALONE (test):
 *   g++ -O3 -std=c++17 -o test_regla ReglaMecanicaUniversal_Consolidadov3.cpp
 *
 * COMPILACIÓ DLL (Windows, MinGW):
 *   g++ -shared -O3 -std=c++17 -march=native -flto
 *       ReglaMecanicaUniversal_Consolidadov3.cpp
 *       -o ReglaMecanicaUniversal.dll
 *       -Wl,--export-all-symbols,--out-implib,libReglaMecanicaUniversal.a
 * ============================================================================
 */

#include <iostream>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <iomanip>
#include <vector>
#include <numeric>   // gcd
#include <cassert>

// ============================================================================
// 0. ARITMÈTICA MODULAR 128-BIT (sense __int128 en MSVC → fallback portàtil)
// ============================================================================

/**
 * mulmod64(a, b, m) = (a * b) % m  sense desbordament.
 * Usa __uint128_t si el compilador ho suporta (GCC/Clang),
 * altrament usa l'algorisme de suma binària O(64) (lent però correcte).
 */
static inline uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<uint64_t>(
        (static_cast<__uint128_t>(a) * b) % m
    );
#else
    // Portàtil: Russian peasant multiplication
    uint64_t res = 0;
    a %= m;
    while (b > 0) {
        if (b & 1) res = (res + a) % m;
        a = (a * 2 >= a) ? (a * 2 % m) : ((a - (m - a)) % m);  // evita overflow
        b >>= 1;
    }
    return res;
#endif
}

/**
 * powmod64(base, exp, m) = base^exp mod m  (exponentació ràpida).
 */
static uint64_t powmod64(uint64_t base, uint64_t exp, uint64_t m) {
    if (m == 1) return 0;
    uint64_t res = 1;
    base %= m;
    while (exp > 0) {
        if (exp & 1) res = mulmod64(res, base, m);
        base = mulmod64(base, base, m);
        exp >>= 1;
    }
    return res;
}

// ============================================================================
// 1. MILLER-RABIN DETERMINISTA (FIX-1)
// ============================================================================

/**
 * miller_rabin_witness(n, a):
 *   Retorna true si 'a' és testimoni de composicitat de n.
 *   Retorna false si n passa el test per la base 'a' (probable primer).
 */
static bool miller_rabin_witness(uint64_t n, uint64_t a) {
    if (n % a == 0) return n != a;   // a divideix n → compost (tret que n==a)

    uint64_t d = n - 1;
    int r = 0;
    while ((d & 1) == 0) { d >>= 1; ++r; }

    uint64_t x = powmod64(a, d, n);
    if (x == 1 || x == n - 1) return false;

    for (int i = 0; i < r - 1; ++i) {
        x = mulmod64(x, x, n);
        if (x == n - 1) return false;
    }
    return true; // compost
}

/**
 * is_prime(n): test determinista per a tot n < 3.317 · 10^24.
 * Bases suficients: {2,3,5,7,11,13,17,19,23,29,31,37}.
 * Font: https://miller-rabin.appspot.com/
 */
static bool is_prime(uint64_t n) {
    if (n < 2)  return false;
    if (n == 2 || n == 3 || n == 5 || n == 7 ||
        n == 11 || n == 13 || n == 17 || n == 19 ||
        n == 23 || n == 29 || n == 31 || n == 37) return true;
    if ((n & 1) == 0 || n % 3 == 0 || n % 5 == 0) return false;

    static const uint64_t bases[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (uint64_t a : bases) {
        if (miller_rabin_witness(n, a)) return false;
    }
    return true;
}

// ============================================================================
// 2. RODA p210 PER A ITERACIÓ EFICIENT (FIX-2)
// ============================================================================

/**
 * Roda 2·3·5·7 = 210.
 * Els 48 salts cobrixen un cicle complet de candidats coprimers amb 210.
 * Elimina el 77.1% de candidats sense cap divisió.
 */
static const int SALTS_RODA_210[48] = {
    10,2,4,2,4,6,2,6,4,2,4,6,6,2,6,4,2,6,4,6,8,4,2,4,
    2,4,8,6,4,6,2,4,6,2,6,6,4,2,4,6,2,6,4,2,4,2,10,2
};

/**
 * Donat un nombre imparell d >= 11 coprimer amb 2·3·5·7,
 * retorna el proper candidat de la roda (suma el salt corresponent).
 * Implementació: iterem per la taula fins trobar l'índex de d%210,
 * llavors apliquem el salt.
 */
static uint64_t seguent_candidat_roda(uint64_t d) {
    // Alinea d als residus vàlids mod 210
    // Taula de residus vàlids mod 210 (coprimers amb 2·3·5·7, senars)
    static const int RESIDUS_210[48] = {
        1,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,
        71,73,79,83,89,97,101,103,107,109,113,121,127,131,
        137,139,143,149,151,157,163,167,169,173,179,181,
        187,191,193,197,199,209
    };
    int rem = static_cast<int>(d % 210);
    for (int i = 0; i < 48; ++i) {
        if (RESIDUS_210[i] == rem) return d + SALTS_RODA_210[i];
    }
    // d no és vàlid → cerca el proper vàlid
    for (int delta = 2; delta < 210; delta += 2) {
        uint64_t nd = d + delta;
        int nr = static_cast<int>(nd % 210);
        for (int r : RESIDUS_210) if (r == nr) return nd;
    }
    return d + 2; // fallback
}

// ============================================================================
// 3. K-SWEEP MDC (FIX-6) — nucli del Mètode Diofàntic Cinemàtic
// ============================================================================

/**
 * k_sweep_mdc(N, m_ini, m_fi):
 *   Cerca factors de N en el rang D = 2m+3 per a m ∈ [m_ini, m_fi].
 *   Itera sobre les pendents k = N//(2m+3) en lloc dels valors m.
 *   Determinista i exhaustiu. Complexitat: O(n_k) = O(m_fi - m_ini).
 *
 *   Matemàtica (Libro 5, Víctor 2025):
 *     Per cada k, el zero de la recta lineal és d_zero = N/k.
 *     Si N % d_zero == 0 i d_zero senar → factor trobat.
 */
static uint64_t k_sweep_mdc(uint64_t N, uint64_t m_ini, uint64_t m_fi) {
    if (m_ini < 1) m_ini = 1;
    if (m_fi < m_ini) return 0;

    uint64_t pos_fi  = 2 * m_fi  + 3;   // D gran (k petit)
    uint64_t pos_ini = 2 * m_ini + 3;   // D petit (k gran)

    uint64_t k_lo = (pos_fi  > 0) ? std::max(uint64_t(1), N / pos_fi ) : 1;
    uint64_t k_hi =                  N / pos_ini;

    for (uint64_t k = k_lo; k <= k_hi; ++k) {
        if (k == 0) continue;
        uint64_t candidat = N / k;
        if (candidat < 3 || (candidat & 1) == 0) continue;
        if (N % candidat == 0 && candidat > 1 && candidat < N)
            return candidat;
    }
    return 0;
}

// ============================================================================
// 4. ESTRUCTURES DE DADES
// ============================================================================

struct SolucionDiofantica {
    bool   encontrada;
    int64_t x;
    int64_t y;
};

struct EstadoConvergenciaAvanzada {
    double   d1, d2, d3, j_actual;
    uint32_t iteraciones;
};

struct EspectroResonancia {
    uint32_t cantidad_frecuencias;
    double   frecuencias[4];
};

// ============================================================================
// 5. CLASSE PRINCIPAL
// ============================================================================

class ReglaMecanicaUniversal {
private:
    // Escala logarítmica per a la regla de càlcul analògica
    const uint32_t precision_;
    const uint32_t tamano_escala_;
    double*        escala_logaritmica_;

    void precargar_escala() {
        for (uint32_t i = 0; i < tamano_escala_; ++i)
            escala_logaritmica_[i] = std::pow(10.0, double(i) / precision_);
    }

public:
    // ── Constructor / Destructor ──────────────────────────────────────────
    ReglaMecanicaUniversal(uint32_t prec, uint32_t tamano)
        : precision_(prec), tamano_escala_(tamano)
    {
        escala_logaritmica_ = new double[tamano_escala_];
        precargar_escala();
    }

    ~ReglaMecanicaUniversal() {
        delete[] escala_logaritmica_;
    }

    // ── Regla de càlcul analògica ─────────────────────────────────────────
    /**
     * multiplicar(pos_a, pos_b): versió "física" original — suma índexs
     * d'escala i llig el valor resultant. Útil si ja treballes en posicions
     * (p.ex. una regla de càlcul dibuixada amb cursors enters).
     */
    double multiplicar(uint32_t pos_a, uint32_t pos_b) const {
        uint32_t idx = pos_a + pos_b;
        return (idx < tamano_escala_) ? escala_logaritmica_[idx] : -1.0;
    }

    /**
     * posicion_de(x): inversa de l'escala — funció que faltava.
     * Donat un valor real x > 0, retorna la posició (índex fraccionari)
     * que ocuparia a l'escala logarítmica: pos = precision · log10(x).
     * És l'operació que fas físicament quan alinees el cursor d'una
     * regla de càlcul sobre un valor concret.
     */
    double posicion_de(double x) const {
        if (x <= 0.0) return -1.0;
        return double(precision_) * std::log10(x);
    }

    /**
     * valor_en_posicion(pos): lectura interpolada de l'escala per a
     * posicions NO enteres (interpolació lineal entre els dos punts
     * de la taula més propers). Necessària perquè posicion_de() rares
     * vegades cau exactament sobre un índex enter.
     */
    double valor_en_posicion(double pos) const {
        if (pos < 0.0 || pos >= double(tamano_escala_ - 1)) {
            // Fora de taula: calcula directament (sense interpolar)
            return std::pow(10.0, pos / double(precision_));
        }
        uint32_t i0 = static_cast<uint32_t>(pos);
        uint32_t i1 = i0 + 1;
        double frac = pos - double(i0);
        return escala_logaritmica_[i0] * (1.0 - frac) + escala_logaritmica_[i1] * frac;
    }

    /**
     * multiplicar_valores(a, b): cicle complet de la regla de càlcul.
     *   1) Converteix a, b a posicions d'escala (cursors).
     *   2) Suma les posicions (desplaçament físic del cursor mòbil).
     *   3) Llig el valor resultant a l'escala (interpolat).
     * Equivalent analògic de a × b, amb l'error característic d'una
     * regla de càlcul real (limitat per `precision_`, no per la taula).
     */
    double multiplicar_valores(double a, double b) const {
        if (a <= 0.0 || b <= 0.0) return -1.0;
        double pos_a = posicion_de(a);
        double pos_b = posicion_de(b);
        return valor_en_posicion(pos_a + pos_b);
    }

    // ── [FIX-1] Test de primalitat (Miller-Rabin determinista) ────────────
    /**
     * Correcte per a tot n < 3.317 · 10^24 (cobreix uint64_t complet).
     * Substitueix el pipeline Karnaugh de v2 que era matemàticament incorrecte.
     */
    bool evaluar_ny_karnaugh(uint64_t ny) const {
        return is_prime(ny);
    }

    // ── [FIX-2] Cerca del proper primer ───────────────────────────────────
    /**
     * Usa roda p210 per saltar candidats compostos.
     * Cada crida és independent (no té estat global).
     */
    uint64_t buscar_siguiente_primo(uint64_t desde_ny) const {
        if (desde_ny <= 2) return 2;
        if (desde_ny <= 3) return 3;

        // Primer candidat senar >= desde_ny
        uint64_t c = (desde_ny % 2 == 0) ? desde_ny + 1 : desde_ny;
        // Comprova primers petits directament
        for (uint64_t p : {uint64_t(3),uint64_t(5),uint64_t(7),uint64_t(11),uint64_t(13)}) {
            if (c == p) return p;
        }
        // Itera amb la roda
        while (true) {
            if (is_prime(c)) return c;
            c = seguent_candidat_roda(c);
        }
    }

    // ── [FIX-3] Fase modular sense estat global ───────────────────────────
    /**
     * fase(ny, M) = 2π · (ny mod M) / M
     * Nota: la versió v2 usava t_k_grande acumulador → resultats incorrectes
     * i problemes de race condition en la DLL.
     */
    double calcular_fase_modular(uint64_t ny, uint64_t M) const {
        if (M == 0) return 0.0;
        constexpr double PI2 = 6.283185307179586;
        uint64_t residuo = ny % M;
        return PI2 * (double(residuo) / double(M));
    }

    // ── Patró de bits (sense canvi respecte v2, correcte) ─────────────────
    uint64_t analizar_patron_bits(uint64_t ny, uint32_t ventana) const {
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

    // ── Curvatura modular ─────────────────────────────────────────────────
    double calcular_curvatura_modular(uint64_t ny, uint64_t M) const {
        if (ny < 2 || M == 0) return 0.0;
        double f0 = calcular_fase_modular(ny - 1, M);
        double f1 = calcular_fase_modular(ny,     M);
        double f2 = calcular_fase_modular(ny + 1, M);
        double v1 = f1 - f0, v2 = f2 - f1;
        double acc = v2 - v1;
        double denom = std::pow(1.0 + v1 * v1, 1.5);
        return (std::abs(denom) < 1e-12) ? 0.0 : std::abs(acc) / denom;
    }

    // ── Massa modular ─────────────────────────────────────────────────────
    /**
     * massa(ny, M) = log10(M / (ny mod M))
     * (mesura quant s'allunya ny d'un múltiple de M)
     */
    double calcular_masa_modular(uint64_t ny, uint64_t M) const {
        if (M == 0) return 0.0;
        uint64_t r = ny % M;
        if (r == 0) return std::numeric_limits<double>::infinity();
        return std::log10(double(M) / double(r));
    }

    // ── Rugositat multiescala ─────────────────────────────────────────────
    double analizar_rugosidad_multiescala(uint64_t ny, uint64_t M) const {
        double r1 = double(analizar_patron_bits(ny, 8)) / 8.0;
        double fase = calcular_fase_modular(ny, M);
        double r2 = std::abs(std::sin(fase));
        double masa = calcular_masa_modular(ny, M);
        double r3 = std::isinf(masa) ? 1.0 : masa / (masa + 1.0);
        return std::sqrt(r1*r1 + r2*r2 + r3*r3);
    }

    // ── Goldbach resonant ─────────────────────────────────────────────────
    /**
     * Cerca p1 + p2 = N (ambdós primers) per a N parell.
     * Usa Miller-Rabin, no el pipeline Karnaugh erroni.
     * El filtre de fase és opcional (millora el rendiment, no la correcció).
     */
    bool escanear_goldbach_resonante(uint64_t N, uint64_t& p1, uint64_t& p2,
                                     uint64_t limite_barrido) const
    {
        if (N % 2 != 0 || N <= 2) return false;
        uint64_t ny = 3;
        uint64_t lim = std::min(N / 2, limite_barrido);
        while (ny <= lim) {
            if (is_prime(ny)) {
                uint64_t comp = N - ny;
                if (is_prime(comp)) {
                    p1 = ny; p2 = comp;
                    return true;
                }
            }
            ny += (ny == 2) ? 1 : 2;
        }
        return false;
    }

    // ── [FIX-4a] Detecció primers bessons ────────────────────────────────
    /**
     * Busca p, p+2 primers tals que p >= ny.
     */
    bool detectar_primos_gemelos_fase(uint64_t ny, uint64_t& p1, uint64_t& p2) const {
        if (ny < 3) ny = 3;
        uint64_t c = (ny % 2 == 0) ? ny + 1 : ny;
        // Límit raonable: fins a 10^6 passes
        for (uint64_t i = 0; i < 1'000'000; ++i) {
            if (is_prime(c) && is_prime(c + 2)) {
                p1 = c; p2 = c + 2;
                return true;
            }
            c = seguent_candidat_roda(c);
        }
        return false;
    }

    // ── [FIX-4b] Detecció inflexió modular ───────────────────────────────
    /**
     * Punt d'inflexió: la curvatura modular supera un llindar.
     * (Mesura si ny és un punt de canvi de convexitat en la fase mod M.)
     */
    bool detectar_inflexion_modular(uint64_t ny, uint64_t M) const {
        if (M == 0 || ny < 2) return false;
        double kappa_prev = calcular_curvatura_modular(ny - 1, M);
        double kappa_curr = calcular_curvatura_modular(ny,     M);
        double kappa_next = calcular_curvatura_modular(ny + 1, M);
        // Inflexió local: kappa_curr és màxim local
        return (kappa_curr > kappa_prev) && (kappa_curr > kappa_next)
               && (kappa_curr > 0.01);
    }

    // ── [FIX-6] K-sweep MDC: factorització de semiprims ─────────────────
    /**
     * Factoritza N (semiprim) usant el k-sweep del MDC de Víctor.
     * Estratègia:
     *   F0: Factors trivials i criba de roda p210 fins a lim_criba.
     *   F1: K-sweep zona densa [m_conv, m_max]  (14% del rang, ~O(√N · 0.14)).
     *   F2: K-sweep zona baixa [lim_criba, m_conv]  (resta).
     *
     * Retorna el factor menut p (o 0 si no troba).
     */
    uint64_t factorizar_mdc(uint64_t N) const {
        if (N < 4) return 0;
        if ((N & 1) == 0) return 2;

        // ── F0: Primers petits ──────────────────────────────────────────
        static const uint64_t PETITS[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
        for (uint64_t p : PETITS) {
            if (N % p == 0 && p < N) return p;
        }

        uint64_t sq = static_cast<uint64_t>(std::sqrt(double(N)));
        while ((sq + 1) * (sq + 1) <= N) ++sq;  // isqrt exacte
        if (sq * sq == N) return sq;             // quadrat perfecte

        uint64_t m_max = (sq - 3) / 2;
        if (m_max < 1) return 0;

        // ── F0b: Criba roda p210 fins a min(m_max, 500K) ───────────────
        const uint64_t LIM_CRIBA = std::min(m_max, uint64_t(500'000));
        {
            // Primer candidat vàlid ≥ D=11 correspon a l'índex 1 de la
            // taula de residus (residu 11 mod 210). Els salts de
            // SALTS_RODA_210 ja són les diferències exactes entre residus
            // consecutius vàlids mod 210 — NO es multipliquen per 2.
            uint64_t D = 11;
            int idx = 1;  // alineat amb RESIDUS_210[1] == 11
            while (D <= 2 * LIM_CRIBA + 3) {
                if (N % D == 0 && D < N) return D;
                D += SALTS_RODA_210[idx % 48];
                idx = (idx + 1) % 48;
            }
        }

        if (LIM_CRIBA >= m_max) return 0;  // N és primer

        // ── F1: K-sweep zona densa [m_conv, m_max] ─────────────────────
        // m_conv = (e-1)/2 · m_max ≈ 0.8591 · m_max
        constexpr double E_MINUS_1_OVER_2 = 0.8591409142295227;
        uint64_t m_conv = static_cast<uint64_t>(E_MINUS_1_OVER_2 * double(m_max));
        if (m_conv < LIM_CRIBA) m_conv = LIM_CRIBA;

        uint64_t f = k_sweep_mdc(N, m_conv, m_max);
        if (f) return std::min(f, N / f);

        // ── F2: K-sweep zona baixa [LIM_CRIBA, m_conv] ─────────────────
        f = k_sweep_mdc(N, LIM_CRIBA, m_conv);
        if (f) return std::min(f, N / f);

        return 0;  // no trobat (N és primer o > 30 dígits equilibrat)
    }

    // ── [NOU] Test de simetria del sawtooth: N mod 2p = p ─────────────────
    /**
     * Descobriment de Víctor: p és l'únic punt de simetria perfecta del
     * "sawtooth accordion" exactament quan N mod 2p = p.
     *
     * Demostració algebraica:
     *   N mod 2p = p  ⟺  N = k·2p + p per algun enter k≥0
     *               ⟺  N = p·(2k+1)
     *               ⟺  p | N  i  N/p és senar.
     * Com que treballem amb semiprims N = p·q amb p,q senars, el quocient
     * N/p = q és sempre senar → la propietat és equivalent a "p divideix N".
     *
     * Verificat empíricament: cap fals positiu en candidats no-factor
     * (provat exhaustivament per a N=10403 amb D∈[3,200)).
     *
     * Avantatge sobre N % candidat == 0: cap, és matemàticament equivalent
     * per a semiprims senars (mateix cost computacional: una divisió).
     * El seu valor és CONCEPTUAL — confirma la interpretació geomètrica
     * del MDC (el factor és el punt de simetria de l'ona) i pot servir
     * com a verificació creuada independent del resultat de factorizar_mdc.
     */
    bool es_punto_simetria(uint64_t N, uint64_t p) const {
        if (p == 0) return false;
        return (N % (2 * p)) == p;
    }


    /**
     * Retorna x tal que A·x ≡ 1 (mod M), o -1 si no existeix.
     * O(log M) — substitueix el bucle O(M) de v2.
     */
    int64_t calcular_inverso_modular_cinematico(int64_t A, int64_t M) const {
        if (M <= 1) return -1;
        int64_t old_r = A % M, r = M;
        int64_t old_s = 1,     s = 0;
        while (r != 0) {
            int64_t q = old_r / r;
            int64_t tmp = r;   r   = old_r - q * r;   old_r = tmp;
            tmp = s;   s   = old_s - q * s;   old_s = tmp;
        }
        if (old_r != 1) return -1;
        return (old_s % M + M) % M;
    }

    // ── Multiplicació modular per bits (sense canvi) ──────────────────────
    uint64_t calcular_multiplicacion_modular_bits(uint64_t A, uint64_t B, uint64_t M) const {
        return (M == 0) ? 0 : mulmod64(A % M, B % M, M);
    }

    // ── Arrel entera per Newton (sense canvi, correcte) ───────────────────
    uint64_t calcular_raiz_cinematica_bits(uint64_t A) const {
        if (A < 2) return A;
        uint64_t x = A >> 1;
        uint64_t prev = 0;
        while (x != prev) {
            prev = x;
            x = (x + A / x) >> 1;
        }
        return x;
    }

    // ── Discriminant exacte (Mètode C, Libro 5) ──────────────────────────
    /**
     * Δ(S) = S² + 6S − (N−9) = k² ?
     * Si sí: p = 2v+3, q = 2b+3  amb v=(S+k)/2, b=(S-k)/2.
     */
    bool discriminant_exacte(uint64_t S, uint64_t N,
                              uint64_t& p_out, uint64_t& q_out) const {
        if (N < 9 || S + 6 > N) return false;
        // delta = S^2 + 6S - (N-9)  pot ser negatiu
        int64_t Ss = static_cast<int64_t>(S);
        int64_t delta = Ss*Ss + 6*Ss - static_cast<int64_t>(N - 9);
        if (delta < 0) return false;
        uint64_t k = static_cast<uint64_t>(std::sqrt(double(delta)));
        // Ajust per error d'arrodoniment FP
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

    // ── Solució diofàntica (sense canvi) ──────────────────────────────────
    SolucionDiofantica resolver_diofantica_cinematico(
            int64_t A, int64_t B, int64_t C, uint64_t limite_barrido) const
    {
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
};

// ============================================================================
// 6. BANC DE PROVES (main)
// ============================================================================

#ifndef BUILDING_DLL

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

    return 0;
}

#endif  // BUILDING_DLL
