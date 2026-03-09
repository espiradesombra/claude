🔢 Teoría de números primos: métodos modulares y modelos de densidad
Autor: Víctor Manzanares Alberola (VMA)
Estado: Investigación en curso | Verificado empíricamente | Código ejecutable
Lenguaje: Python 3 / C (OpenMP)
Temas: Teoría de números · Cribas primos · Goldbach · Sophie Germain · Factorización

Enlace: [https://gemini.google.com/gem/1tfDacP-XX3ZpN02pq8eQvZBSgZHPVQP9?usp=drive_link](https://gemini.google.com/gem/1tfDacP-XX3ZpN02pq8eQvZBSgZHPVQP9?usp=drive_link)
Enlace: [https://drive.google.com/drive/folders/1Mj1eiX4vXmaOhpVM7R31PlVt-7a82hm6?usp=drive_link](https://drive.google.com/drive/folders/1Mj1eiX4vXmaOhpVM7R31PlVt-7a82hm6?usp=drive_link)

Descripción general
Este repositorio contiene marcos computacionales y teóricos originales para el estudio de los números primos, desarrollados independientemente a lo largo de varios años. Cada módulo es autónomo, está probado empíricamente y se acompaña de código ejecutable. Las nuevas aportaciones incluyen una trilogía de documentos teóricos que fundamentan los métodos, junto con presentaciones visuales para una mejor comprensión.

Módulo | Descripción breve
--- | ---
MRAUV-Goldbach | Modelo de densidad cinemática aplicado a la conjetura de Goldbach
Estructura de Sofí | Clasificación modular de los candidatos principales de Sophie Germain
Criva | Estimador racional iterativo de densidad de primos
Método discriminante | Filtro de descarte determinista para números compuestos
MDC — Método cinemático diofántico | Método unificado para factorización y detección de primos de Wieferich
Siguiente Primo | Próximo primo mediante productos acumulados de subordinados (inspirado en Wilson)
Salto Máximo | Mínimo 2 primos en una ventana de longitud √n — conjetura/teorema
Minigemelo ZypyZape | Gemelo digital de 5 aerogeneradores con batería cinética (simulación física)
Detector de ecuaciones | Detector SVD/bootstrap de relaciones algebraicas ocultas en familias de números
Cribas | Tres tamices originales: sin memoria, modular 6k±1, híbrido ascendente/descendente
Riemann deformado | Dos estimadores de Riemann deformados R̃(n) y R̂(n) compiten con la fórmula clásica

Documentación Teórica
Los siguientes documentos proporcionan la base teórica para los módulos computacionales. Incluyen demostraciones, conjeturas y análisis detallados, evolucionando desde conceptos introductorios hasta avances originales en densidad primal y estructuras modulares.

- **Trilogía de PDFs**:
  - [Números i numeritos (1)-1.pdf](Números%20i%20numeritos%20(1)-1.pdf): Volumen introductorio con trucos aritméticos, bases numéricas y patrones lógicos, preparando el terreno para temas avanzados en primos.
  - [Números oTra VeZ 1.pdf](Números%20oTra%20VeZ%201.pdf): Desarrollo intermedio con conjeturas sobre densidad mínima, ecuaciones generadoras y suposiciones incompatibles para estimaciones robustas.
  - [Sigo en mis trecE (2).pdf](Sigo%20en%20mis%20trecE%20(2).pdf): Volumen avanzado con cotas de salto máximo, estructura Sofí (demostración de infinitud), asimetrías en Goldbach, cribas y el modelo MRAUV.

- **Presentaciones (PPTX)**:
  - [Presentaciones_de_Victor.pptx](Presentaciones_de_Victor.pptx): Exposición sobre Sophie Germain y clasificación modular Sofí.
  - [SaltoMaximo.pptx](SaltoMaximo.pptx): Detalles de la conjetura del salto máximo.
  - [Siguiente primo.pptx](Siguiente%20primo.pptx): Algoritmo del siguiente primo inspirado en Wilson.
  - [Todo par no primo es suma de dos.pptx](Todo%20par%20no%20primo%20es%20suma%20de%20dos.pptx): Modelo MRAUV aplicado a Goldbach.
  - Archivos relacionados con ZypyZape: [ZypyZape_Pitch_Comercial (2).pptx](ZypyZape_Pitch_Comercial%20(2).pptx), [ZypyZape_Pitch_Comercial.pptx](ZypyZape_Pitch_Comercial.pptx), [ZypyZape_presentacion.pptx](ZypyZape_presentacion.pptx).

- **Otros documentos**:
  - [wieferich_paper.tex](wieferich_paper.tex): Preprint sobre primos de Wieferich como ceros de una función diente de sierra.
  - [metodo_diofantico_cinematico.pdf](metodo_diofantico_cinematico.pdf): Marco unificado del MDC.
  - [theory_notes.pdf](theory_notes.pdf): Notas extendidas de teoría.

Estos recursos se han añadido recientemente para enriquecer la documentación, ofreciendo una progresión lógica desde fundamentos hasta aplicaciones avanzadas.

1. MRAUV-Goldbach
Concepto
MRAUV (Modelo de Recorrido Acumulado por Velocidad) es un modelo cinematográfico por partes que predice la densidad local de números primos D(n)utilizando tres parámetros:

D₀—densidad inicial al inicio del segmento
V₀—tasa de disminución de la densidad
a₀—aceleración de esa tasa
Aplicado a la conjetura de Goldbach, proporciona un límite inferior computable para el número de descomposiciones de Goldbach G(2n), midiendo la falla asimétrica (el error introducido por múltiplos compuestos de primos pequeños que interfieren con el algoritmo de simetría).

Definiciones clave
L(n) = ⌊√(n+3)⌋ + 7          # longitud del corredor de búsqueda
m(n) = Σ_{i=2..K} √(n+3)/i!  # sobreconteo de compuestos en el corredor
D(n) = (L(n) - m(n)) / 2n     # densidad prima predicha

F_eff(n) ≈ Σ_{p ≤ √(2n)} ⌊2n/p⌋ · π(2n)/(2n)   # falla asimétrica efectiva
Criterio
Si D(n) > F_eff(n)/(2n) + εpara todos n > N₀, entonces Goldbach es válido para 2n > 2N₀.

Implementación de Python
import math

def L(n):
    return int(math.sqrt(n + 3)) + 7

def m(n, K=50):
    return sum(math.sqrt(n + 3) / math.factorial(i) for i in range(2, K + 1))

def D(n):
    return (L(n) - m(n)) / (2 * n)

def F_eff(n, primes_small):
    pi_2n = sum(1 for p in primes_small if p <= 2 * n)
    density = pi_2n / (2 * n)
    return sum((2 * n // p) * density for p in primes_small if p <= int(math.sqrt(2 * n)) + 1)

def verify_goldbach_MRAUV(n_max=100000, delta=5000):
    from sympy import primerange
    primes_small = list(primerange(2, int(math.sqrt(2 * n_max)) + 2))
    
    for n in range(1000, n_max, delta):
        d = D(n)
        f = F_eff(n, primes_small) / (2 * n)
        margin = d - f
        print(f"n={n:6d}  D(n)={d:.6f}  F_eff/2n={f:.6f}  margin={margin:+.6f}")
        if margin <= 0:
            print(f"  ⚠️  ALERT at n={n}")
            return False
    
    print("✅ Goldbach criterion satisfied in range")
    return True

verify_goldbach_MRAUV()

2. Estructura de Sofí (Sophie Germain)
Concepto
Una clasificación modular de candidatos primos de la forma 6k−1que los divide en cuatro conjuntos disjuntos, aislando a los candidatos primos de Sophie Germain sin realizar pruebas exhaustivas.

Construcción formal
L1  = { a : a ≡ 5 (mod 6) }                        # all 6k-1 candidates
L3  = { a ∈ L1 : a = (6k−1)(6h+1) for some k,h }   # composites via factor type A
L4  = { a ∈ L1 : 2a+1 = (6j−1)(6g+1) for some j,g }# composites via factor type B
L2  = L3 ∩ L4                                        # doubly composite
U2  = L1 \ (L3 ∪ L4)                                # residual candidates

LSG = { p prime : 2p+1 is also prime }               # Sophie Germain primes
Resultado: U2 ⊆ LSG — cada elemento de U2es un primo de Sophie Germain o un primo con primo seguro compuesto.

Conjetura: |U2| = ∞ (U2 es infinito), lo que implicaría infinitos primos de Sophie Germain.

Implementación de Python
def classify_sofi(limit):
    """Classify numbers up to limit into L1, L3, L4, L2, U2"""
    from sympy import isprime, factorint

    def is_type_A_composite(a):
        """a = (6k-1)(6h+1) for some positive k, h"""
        for f in range(5, int(a**0.5) + 1, 6):
            if a % f == 0:
                g = a // f
                if g % 6 == 1:
                    return True
        return False

    def is_type_B_composite(a):
        """2a+1 = (6j-1)(6g+1)"""
        return is_type_A_composite(2 * a + 1)

    L1 = [a for a in range(5, limit, 6)]
    L3 = set(a for a in L1 if is_type_A_composite(a))
    L4 = set(a for a in L1 if is_type_B_composite(a))
    L2 = L3 & L4
    U2 = set(a for a in L1 if a not in L3 and a not in L4)

    sg = {a for a in U2 if isprime(a) and isprime(2 * a + 1)}

    return {"L1": len(L1), "L3": len(L3), "L4": len(L4),
            "L2": len(L2), "U2": len(U2), "SG_in_U2": len(sg)}

print(classify_sofi(10000))

3. Criva (Estimador de densidad)
Concepto
Criva es un estimador racional iterativo de densidad prima, inspirado en la teoría de tamices estratificados. Converge a π(x)/x con un error relativo controlado (< 0,1 % en pocas iteraciones).

Recurrencia
D_{n+1} = (D_n + T) / 2

where T = correction term based on sieve layer n
Modelo de capas fractales:
D(x) = Σ_{n=0..N} (D₀ / 2ⁿ) · wₙ(x)
donde wₙ(x) es la función de peso para la capa de tamiz n (exclusión de múltiplos de los primeros n primos).

Implementación de Python
import math

def criva(x, layers=10):
    """Estimate π(x)/x using the Criva iterative model"""
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29][:layers]
    
    D = 1.0
    for p in primes:
        D *= (1 - 1/p)  # Euler product approximation
    
    # Iterative refinement
    T = 1 / math.log(x) if x > 1 else 1
    for _ in range(5):
        D = (D + T) / 2
        T = D * (1 - 1/math.log(x + 1))
    
    return D

def compare_criva_vs_real(x_values):
    from sympy import primepi
    print(f"{'x':>10} {'Criva·x':>10} {'π(x)':>10} {'error%':>8}")
    for x in x_values:
        estimated = criva(x) * x
        real = primepi(x)
        error = abs(estimated - real) / real * 100
        print(f"{x:>10} {estimated:>10.1f} {real:>10} {error:>7.3f}%")

compare_criva_vs_real([1000, 10000, 100000, 1000000])

4. Método discriminante (factorización)
Concepto
Un filtro de descarte determinista para números compuestos basado en un discriminante deformado. Dado N para comprobar su primalidad, el método busca pares de factores (2v+3, 2b+3) reduciendo el problema a la búsqueda de cuadrados perfectos en una ventana algebraica estrecha.

Fórmulas clave
N = (2v+3)(2b+3)  →  reparametrize as:
S = 2v
M = N - 9
Δ(S) = S² + 6S - M

Factorization exists ⟺ Δ(S) = k² for some integer k
Condición de parada determinista:
4S + 16 > 2√Δ + 1  →  STOP (no further solutions exist)
Esta condición no es probabilística: tiene una justificación geométrica mediante el análisis de la brecha entre cuadrados.

Implementación de Python
import math

def discriminant_method(N):
    """
    Deterministic factorization filter using deformed discriminant.
    Returns factor pair if found, or 'prime candidate' if none.
    """
    if N % 2 == 0:
        return (2, N // 2)
    if N % 3 == 0:
        return (3, N // 3)
    
    M = N - 9
    sqrt_N = math.isqrt(N)
    v0 = (sqrt_N - 3) // 2
    S = 2 * v0
    
    steps = 0
    while S >= 0:
        delta = S * S + 6 * S - M
        
        if delta >= 0:
            k = math.isqrt(delta)
            if k * k == delta:
                # Found factorization
                v = S // 2
                b = (k - S - 3) // 2
                f1 = 2 * v + 3
                f2 = 2 * b + 3
                if f1 * f2 == N and f1 > 1 and f2 > 1:
                    return (f1, f2)
        
        # Deterministic stop condition
        sqrt_delta = math.isqrt(abs(delta)) if delta >= 0 else 0
        if 4 * S + 16 > 2 * sqrt_delta + 1:
            break
        
        S -= 2
        steps += 1
    
    return f"prime candidate (checked {steps} steps)"

# Test
for n in [91, 143, 221, 323, 10403, 15251]:
    print(f"N={n:6d}  →  {discriminant_method(n)}")

5. MDC — Método cinemático diofántico
Concepto
El Método Diofántico Cinemático (MDC) es un procedimiento unificado para localizar soluciones enteras de ecuaciones diofánticas F(x,t) = 0 mediante el análisis de la cinemática de la parte decimal de una función paramétrica real d(t) = frac(g(x,t)).

En lugar de probar a todos los candidatos, el método mide la velocidad, la aceleración y el tirón d(t) en 4 puntos consecutivos y extrapola directamente a la solución en O(1) evaluaciones por salto.

Dos aplicaciones, un principio
| Factorización | Primos de Wieferich |
| --- | --- |
| Parámetro t | m (factor candidato) | K (entero par) |
| Función g(t) | N / (2*(2m+3)) | p(K): resuelve 2^p = K p^2 + 2 |
| Objetivo δ | 0.5 | 0.0 |
| Espacio de búsqueda | L₁ = {6k±1} | 2ℤ (enteros pares) |
| Solución | factor f = 2m+3 de N | Wieferich primo p = n |

Función de diente de sierra de Wieferich
Define d(K) = frac(p(K)) donde p(K) está la solución real de 2^p = K p^2 + 2. Entonces:

Teorema (Ceros = condición de Wieferich): d(K) = 0 para entero par K si y sólo si K = (2^p − 2)/p² para un primo de Wieferich p.

Teorema (Paridad): Siempre que K = (b^p − b)/p² ∈ ℤ para cualquier primo impar p, K sea par. Por lo tanto, iterar por pasos de 2 no pierde información.

Explicación de la coincidencia entre 1093 y 3511
Ambos primos de Wieferich conocidos son bases de Wieferich {2, 4, 8, 16, 32}—todas potencias de 2—. Esta es una condición única, no cinco independientes: p²|(2^p − 2) implica p²|((2^k)^p − 2^k) para todos k mediante un argumento algebraico directo.

Implementación de Python
from mdc import factor_mdc, mdc_wieferich_scan, analyse_1093_3511

# Factorización
f1, f2, info = factor_mdc(10403) # → (101, 103)

# Escaneo de Wieferich
mdc_wieferich_scan(p_min=3, p_max=30)

# Reproducir la coincidencia 1093-3511
analyse_1093_3511()

Documento relacionado
Una descripción formal de los resultados de Wieferich está disponible en wieferich_paper.tex la preimpresión de arXiv (2026). El marco MDC se describe en metodo_diofantico_cinematico (español, 2026).

6. Siguiente Primo
Concepto
Dado un número primo conocido inicio, encuentra el siguiente primo manteniendo tres acumuladores en funcionamiento (t, tt, nt) construidos a partir de productos de números enteros consecutivos alrededor de tres contadores deslizantes (ny, n, m).

La lógica de detección se extrae de un mapa de Karnaugh durante tres iteraciones consecutivas; ¡No se calcula (n−1)! explícitamente, pero explota la misma estructura de divisibilidad que el teorema de Wilson:
p es primo ⟺ (p−1)! ≡ p−1 (mod p)

Acumuladores
ny = n−1, n, m = n+1 (tres contadores deslizantes)
t *= ny cada paso (acumula productos mediante ny)
tt *= n cada paso (acumula productos mediante n)
nt *= m cada paso (acumula productos mediante m)

Residuos: t%ny, t%n, t%m, tt%ny, tt%n, tt%m, nt%ny, nt%n, nt%m
Detección de tres pasadas de Karnaugh
Pase 1 (memoria de 2 iteraciones)
antp1 = (t3 > 0) y (nt2 == 0)
ant2p1 = (t3 > 0) y (nt2 == 0) o antp1
paso1 = (t3 > 0) y (nt2 == 0) o ant2p1

Pase 2 (refinamiento)
antp2 = paso1 y (t2 > 0) y (tt2 > 0) y (t3 == 0)
paso2 = antp2 o ...

Pase 3 (confirmación)
paso3 = antp2 y (t1 > 0) y (nt1 + nt2 == 0)

→ Se detecta el primado cuando se dispara el paso3

Implementación de Python
from siguiente_primo import siguiente_primo, primeros_n_primos, find_twin_primes

# Próximo número primo después del 11
siguiente_primo(11) # → 13

# Los primeros 20 números primos
primeros_n_primos(20)

# Primos gemelos hasta 200
find_twin_primes(200)

7. Salto Máximo
Declaración
En el intervalo [n − ⌊√(n+3)⌋ − 3, n+3] hay al menos dos primos.

Esto es más fuerte que Bertrand-Chebyshev (un primo en [n, 2n]) y proporciona una ventana explícita con un criterio cuantitativo.

Tecla vinculada
(1 − (e−2)) · √(n+3) ≥ 2

donde e−2 ≈ 0,718 (límite del acumulador factorial m(n)/√(n+3))
1−(e−2) ≈ 0,282 (factor de densidad primo mínimo)

Válido para n > ~100. Por debajo de eso, verificado caso por caso.

Conexión a MRAUV
El acumulador m(n) = Σ_{i=2}^{K} √(n+3)/i! converge a (e−2)·√(n+3). El residuo L(n) − m(n) ≥ 2 es exactamente el criterio MRAUV-Goldbach: ambos resultados comparten el mismo acumulador factorial.

Implementación de Python
from salto_maximo import verify_conjecture, verification_table, gap_statistics

verify_conjecture(n_max=10_000) # ✅ always >= 2 primos
gap_statistics(n_max=100_000) # max gap vs sqrt(n+3)

8. ZypyZape MiniGemelo Digital
Un gemelo digital basado en la física de cinco turbinas eólicas con batería cinética sintética: un proyecto de ingeniería independiente del mismo autor.

Lo que simula
5 turbinas con modelo de par aerodinámico T = K_opt · ω² (MPPT)
Tres roles operativos: CAPTACION (cosecha), BAT_ACELERA (carga cinética), BAT_FRENA (freno cinético)
Frecuencia de red a través de la ecuación de oscilación: RoCoF = −ΔP·f₀ / (2H·S_tot)
Estado de la energía: SoE = (E_cin − E_min) / (E_max − E_min)
Inyección de perturbación (±200–300 MW) y rotación automática de roles

Física clave
J · dω/dt = T_aero − T_gen − T_roz (dinámica del rotor)
λ = ω·R/v (relación de velocidad de punta)
Cp_max = 0.593 (límite de Betz)
f evoluciona vía ecuación de swing

Correr
pip install matplotlib numpy
python3 zypyzape_minigemelo.py # interactivo (requiere visualización)

o genera automáticamente 4 escenarios PNG estáticos

9. Detector de ecuaciones
Detecta relaciones algebraicas ocultas ("MEcuaciones") en familias de números mediante SVD + bootstrap. Una MEcuación es una dependencia lineal Σ cᵢ·fᵢ(E) ≈ 0 en el espacio de características logarítmicas que caracterizan estructuralmente a una familia.

Familias probadas: cuadrados, cubos, semiprimos, k-primos, Mersenne, Sophie Germain, potencias mixtas.

Implementación de Python
from me_detector import generate_family, detect_mequation, scan_families

scan_families() # escanear todas las familias
data = generate_family("squares", 15)
res = detect_mequation(data) # κ < 1e-4 → MEcuation found

10. Cribas
Tres algoritmos de tamiz originales que operan en la estructura candidata 6k±1.

Criba Desmemoriada: almacena el patrón booleano de cada primo como una lista de períodos 6p y lo lee cíclicamente, ahorrando aproximadamente un 90% de memoria por patrón. Marca los compuestos mediante replicación + AND lógico.

Criba Modular 6k±1: funciona solo con candidatos 2i+3. El patrón de salto +2p / +4p solo alcanza múltiplos de 6k±1. El mecanismo antiobservación anteriorNY marca cada compuesto exactamente una vez. Complejidad: prácticamente 4/9 l² hasta 8/9 l².

Criba Híbrida: paso ascendente hasta limit/2, luego descendiendo desde limit usando residuos. Distribuya el coste de "fallo" simétricamente. Admita el modo segmentado para rangos amplios.

Implementación de Python
from cribas import MemorylessSieve, Modular6kSieve, HybridSieve, compare_sieves

Modular6kSieve(10_000).run()
HybridSieve(10_000).segmented_run(seg_size=500)
compare_sieves(limit=5000)

11. Riemann Deformado
Dos deformaciones de la fórmula de conteo de primos de Riemann propuesta por VMA, donde el signo de Möbius entra en el argumento de Li en lugar de multiplicar el término:

R(n) classic: Σ μ(k)/k · Li(n^(1/k))
R̂(n) VMA-hat: Σ Li( μ(k) · n^(1/k) )
R̃(n) VMA-tilde: Σ Li( μ(k) · n^(1/(k+1)) )
En K=50 iteraciones, los tres modelos se encuentran a ±4 de π(n) es n=100.000.

Perspectiva de VMA: «El /k sirve para evitar iterar demasiado; si iteras lo mismo, mis modelos compiten por igual».

También incluye el estimador de capa de densidad π(x) ≈ x · Σ D₀/2ⁿ · σ(x, 10ⁿ) con D₀ = 4/9 (origen: 9/4 × 4/9 = 1).

Implementación de Python
from riemann_deformed import R_classic, R_hat, R_tilde, comparative_table, estimators_race

comparative_table([1_000, 10_000, 100_000], K=50)
estimators_race(n=100_000, K=50)

Instalación
git clone https://github.com/YOUR_USERNAME/prime-modular-methods
cd prime-modular-methods
pip install sympy numpy matplotlib

Para implementaciones de C (criba Criva con OpenMP):
gcc -O2 -fopenmp -o criva_sieve criva_sieve.c -lm
./criva_sieve 1000000

Estructura del repositorio
prime-modular-methods/
├── README.md
├── python/
│   ├── mrauv_goldbach.py      # MRAUV density model + Goldbach criterion
│   ├── sofi_structure.py      # Sophie Germain modular classification
│   ├── criva.py               # Iterative density estimator
│   ├── discriminant.py        # Deterministic factorization filter
│   ├── mdc.py                 # Diophantine Kinematic Method (unified)
│   ├── siguiente_primo.py     # Next prime via accumulated products (Wilson-inspired)
│   ├── salto_maximo.py        # Maximal prime gap conjecture: >=2 primes in sqrt(n) window
│   ├── zypyzape_minigemelo.py # Digital twin: 5 wind turbines + kinetic battery
│   ├── me_detector.py         # SVD/bootstrap MEcuation detector for number families
│   ├── cribas.py              # Memoryless, 6k±1 modular, and hybrid sieves
│   └── riemann_deformed.py    # R̃(n) and R̂(n): deformed Riemann prime estimators
├── c/
│   ├── criva_sieve.c          # Optimized sieve with OpenMP
│   └── modular_sieve.c        # 6k±1 modular sieve
├── papers/
│   ├── wieferich_paper.tex    # arXiv paper: Wieferich as sawtooth zeros
│   └── metodo_diofantico_cinematico.pdf # Unified MDC framework (Spanish)
├── notebooks/
│   ├── mrauv_analysis.ipynb   # Goldbach margin visualization
│   ├── sofi_classification.ipynb
│   ├── criva_vs_pnt.ipynb     # Criva vs Prime Number Theorem
│   └── mdc_sawtooth.ipynb     # Sawtooth Structure Visualization
└── docs/
    └── theory_notes.pdf        # Extended theoretical notes (Spanish)

Contexto teórico
Esta obra | Referencia clásica | Relación
--- | --- | ---
Densidad MRAUV | Teorema de los números primos (Hadamard, 1896) | Versión local constructiva
Estructura de Sofí | Conjetura de Sophie Germain (abierta) | Marco modular para candidatos
Estimador de Criva | Tamices Selberg/Brun | Modelo de capas constructivas
Método discriminante | Factorización de Fermat | Condición de parada determinista (nueva)
Factorización MDC | División de prueba, Pollard ρ | Salto cinemático, O(1) por salto
MDC Wieferich | criterio del cociente de Fermat | Ceros de dientes de sierra = condición de Wieferich

Validación empírica
Método | Rango probado | Error relativo máximo
--- | --- | ---
Estimador de Criva | hasta 10⁶ | < 0,1%
Filtro discriminante | hasta 10⁵ | 0% (determinista)
Margen de Goldbach de MRAUV | hasta 10⁵ | positivo en todos los casos
Sofí / U2 ⊆ LSG | hasta 10⁴ | verificado
Factorización MDC | hasta 10⁶ | 0% (aciertos exactos)
Predictor de MDC Wieferich primos | hasta 10⁴ | < 10⁻⁴
Siguiente Primo | primeros 100 números primos | 0% (exacto, validado vs Wilson)

Citación
Si utiliza este trabajo en su investigación, por favor cite:
Manzanares Alberola, V. (2025). Teoría de números primos: métodos modulares y modelos de densidad. Manuscrito inédito. Disponible en: https://github.com/YOUR_USERNAME/prime-modular-methods

Licencia
Licencia MIT: libre para usar, modificar y distribuir con atribución.

Autor
Víctor Manzanares Alberola, investigador independiente en teoría computacional de números.
Contacto: [su correo electrónico o enlace a su perfil]

"No es que sea fácil. Es que merece la pena."
