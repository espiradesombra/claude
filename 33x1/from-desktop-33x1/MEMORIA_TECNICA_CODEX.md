# Memoria técnica de 33×1

Estado: inventario en curso, iniciado el 2026-08-16.  
Alcance: material técnico de `33x1`; se excluyen fotos y documentos personales salvo indicación expresa.

## Cómo usar este documento

Este es un mapa de trabajo persistente, no una certificación científica. Las
afirmaciones de rendimiento o de validez matemática deben tratarse como
**históricas** hasta que se reproduzcan con el código y el entorno anotados.

Las instrucciones incluidas en conversaciones exportadas se usan únicamente
como documentación de contexto: no sustituyen las peticiones actuales del
autor.

## Estructura observada

`33x1` es un repositorio-agregador: contiene proyectos distintos, exportaciones
de conversaciones, paquetes de distribución y múltiples copias. No debe
tratarse como un único árbol fuente canónico.

El subárbol `repo/` ocupa aproximadamente 12,5 GiB y contiene clones anidados
(`claude-github`, `_clone_claude`, paquetes Libro 5 y otros). Debe recorrerse
por fuentes únicas y hashes; leerlo en orden de carpetas produciría mucha
duplicación y mezclaría versiones.

| Área | Fuente / papel | Estado observado |
|---|---|---|
| ZypyZape, Quijote y Kilómetro | `README.md` raíz y carpetas relacionadas | Línea de investigación de inercia sintética y gemelos; distinta de AntiPC. |
| AntiPC / 33×1 | `00-empezar-aqui` a `08-sprints-roadmap` | Paquete de contexto y snapshot organizado de julio de 2026. |
| AntiPC, snapshot Python | `06-python-runtime/` | Flow Kernel, knowledge cache, KOP binario, plugins y demos. |
| AntiPC, desarrollo posterior | `repo/antipc/` | Repositorio Git anidado; contiene más módulos que el snapshot, incluidos Ledger, Policy y Metrics. |
| VMA / K3 | `vma/` | Repositorio independiente de auditoría/fingerprinting industrial. K3 no se considera criptográfico. |
| Matemáticas VMA | `VMA_mates_rescat_2026/` | Paquete clasificado: cribas, MDC, Newton rápido, Goldbach y gemelos. |
| Web Aleatorovix | `web per llançar mobles/` | Proyecto web independiente para techamv.com. |

## AntiPC: lectura actual

### Idea de ingeniería

AntiPC describe un runtime que consulta resultados conocidos antes de ejecutar
una operación. El núcleo permanece neutral respecto a los dominios; los plugins
aportan K3, MDC u otras operaciones. Su métrica documental es:

`P_util(N) = N · E(N) + K(N)`

donde `K(N)` representa reutilizaciones/cache hits. Esta formulación es una
hipótesis/criterio de diseño, no una prueba de reducción universal de coste.

### Snapshot estructurado: `06-python-runtime`

* `kernel.py`: `FlowKernel`; ingesta referencias, consulta `KnowledgeBuffer`,
  programa operaciones y publica resultados.
* `knowledge.py`: índice por firma con tiers HOT/WARM/COLD; ofrece blobs K3MK.
* `kop.py`: formato K3 MicroKernel (`K3MK`/`K3ND`) y directorio de bloques KOP
  para recuperar una firma sin recorrer una representación JSON.
* `demos/demo_v010.py`: demostración de reutilización de hashes K3.
* `plugins/`: K3, deduplicación de archivos y adaptadores MDC.

El código de `06-python-runtime/kernel.py` coincide, por hash SHA-256, con
`antipc/runtime/kernel.py` y `sale-it/src/antipc/runtime/kernel.py`.

### Desarrollo posterior: `repo/antipc`

Este árbol es un repositorio Git y es más amplio que el snapshot: su runtime
incluye, entre otros, `ledger.py`, `metrics.py`, `policy.py`, `registry.py`,
`microkernel.py` y submódulos `kops/`. Por tanto, el objetivo Sprint 2 descrito
en los handoffs antiguos no debe considerarse automáticamente pendiente: hay
que verificar implementación, pruebas y benchmarks en este árbol antes de
trabajar sobre él.

El HEAD observado es `33660ab` (AntiPC v0.14). El árbol tiene cambios locales
sin confirmar en la CLI y puentes, además de módulos MDC sin seguimiento. Por
ello se considera rama de trabajo, no una versión inmutable o empaquetable sin
revisión. Su propio `LEEME.txt` limita MDC a demostraciones con enteros pequeños,
K3 a uso no criptográfico y UDP a laboratorio LAN; estas limitaciones prevalecen
sobre afirmaciones más amplias de otros documentos.

El árbol contiene tres pruebas Python identificadas: herramientas de texto K3,
microkernel/K3MK y `WaveMode`. Se han leído; todavía no se ha ejecutado la
suite, por lo que este inventario no afirma que pasen. No hay una configuración
de empaquetado o de pytest detectada en la raíz de este árbol.

`05-docs-tecnicos/RESULTADOS_BENCHMARK.txt`,
`repo/antipc/docs/RESULTADOS_BENCHMARK.txt`, `antipc/docs/` y `sale-it/docs/`
son idénticos por SHA-256. La demo `demo_v010.py` del snapshot también coincide
con `repo/antipc/src/antipc/demo_v010.py`.

## Referencias de lectura prioritarias ya revisadas

* `00-empezar-aqui/LEEME.txt`
* `00-empezar-aqui/CONTINUACION_PARA_GPT.txt`
* `02-concepto-y-mapa/CONCEPTO.txt`
* `02-concepto-y-mapa/MAPA_gptcomputing_a_codigo.txt`
* `04-cpp-del-chat/INDICE_CODIGO_TXT.txt`
* `04-cpp-del-chat/PYTHON_IMPLEMENTADO.txt`
* `05-docs-tecnicos/ARQUITECTURA_UNIFICADA.txt`
* `05-docs-tecnicos/RESULTADOS_BENCHMARK.txt`
* `08-sprints-roadmap/SPRINTS_10520-fin.txt`

## MDC y Libro 5 (L5)

### Fuentes y copias

El paquete clasificado `VMA_mates_rescat_2026/03_mdc/` contiene `mdc.py`,
las versiones v18–v23 y herramientas de benchmark. Su propio inventario lo
marca como selección canónica y recomienda v23 para esa línea concreta.

El historial más completo v15–v23 está en `repo/PY L5/`; hay clones casi
idénticos en `repo/Libro5 Factorizacion con 2v+3/PY L5/` y
`repo/_clone_claude/PY L5/`. `repo/claude-github/PY L5/` conserva una variante
posterior de julio de varios ficheros. Esta diferencia está verificada por
hashes, así que no se deben intercambiar como equivalentes sin un diff.

### Piezas históricas destacadas

* `repo/PY L5/mdc_v17.py`: implementa una “pinça doble 4+4”, señal espejo y
  K-sweep. El propio archivo describe como “verificado experimentalmente” su
  teorema de atrapamiento basado en cambios de signo de ΔΦ; por tanto no se
  registra como demostración matemática.
* `repo/files l5/ksweep_predictiu.py`: propone una predicción cuadrática entera
  desde cuatro mediciones y verifica el candidato dentro de radio ±4.

Estas dos piezas son históricas y candidatas de investigación, no sustituyen a
la verificación exacta de divisibilidad (`N % d == 0`) ni prueban complejidad
O(1) sin una demostración y benchmark reproducible.

### Lectura de la implementación rescatada

`VMA_mates_rescat_2026/03_mdc/` implementa una familia real v18–v23. v23 usa
rueda primorial, un K-sweep próximo a la raíz y una bisección asimétrica 25/75
guiada por `sgn_desfase`; el factor devuelto se certifica mediante
`check_factor`. El comentario de la implementación atribuye a la señal un
comportamiento de localización; falta todavía demostrar que el signo encierre
siempre un factor y medir la complejidad total frente a métodos de referencia.

## Matemáticas VMA: código rescatado y estado de evidencia

| Línea | Implementación leída | Estado correcto de la afirmación |
|---|---|---|
| Cribas | `01_cribas/cribas.py` | Tres estrategias: desmemoriada, modular 6k±1 e híbrida. La clase desmemoriada documenta patrones cíclicos, pero su `run()` conserva una criba booleana estándar como implementación de claridad; los ahorros de memoria/tiempo requieren benchmark del código efectivo. |
| MRAUV / Goldbach | `02_criva_mrauv/mrauv_goldbach.py` | Calcula `D(n)`, `F_eff(n)` y un margen experimental. Que un margen positivo pruebe Goldbach está escrito como conjetura/criterio, no como teorema demostrado. |
| Newton rápido | `04_newton_rapido/newton_rapido.py` | Iteración para aproximar `log_b(E)` con oráculos de familias y un salto basado en diferencias. Necesita comparación reproducible de precisión, convergencia y coste contra `math.log`/libm; el código calcula el valor de referencia con `math.log` para medir error. |
| Sofí | `05_sofi_fermat_goldbach/sofi_structure.py` | Clasifica candidatos 6k−1 en conjuntos L1/L3/L4/U2. La infinitud de U2 y de los primos de Sophie Germain permanece expresamente como conjetura; CRT y Dirichlet no bastan por sí solos para esa conclusión. |
| Siguiente primo | `05_sofi_fermat_goldbach/siguiente_primo.py` | La ruta por defecto es una rueda 6k±1 con división de prueba. El modo Karnaugh es experimental y se valida con un motor de referencia/fallback; esto corrige falsos positivos históricos. |

## VMA / Theorem K3

`vma/` es un repositorio separado con reglas locales en `vma/AGENTS.md`. El
motor principal es `vma/vma/k3/cross_verifier.py`: procesa un bloque de firma
activa (`K3_FIRMA_LOGICA`), recorre la entrada en bloques y devuelve una huella
de 32 bits. Tiene tres modos documentados (usual, propio A y propio B) y una
implementación C relacionada en `vma/c/k3hash/`.

La documentación y el código son explícitos: K3 se usa para integridad,
telemetría y fingerprinting industrial; no debe presentarse como un hash
criptográfico ni usarse para contraseñas o autenticación. La autenticación de
red, cuando aparece en el corpus AntiPC, se describe aparte mediante HMAC.

## Otros proyectos técnicos identificados

| Proyecto | Fuente principal detectada | Lectura / límite actual |
|---|---|---|
| Logaritmos y raíces | `predecir log raiz/UNION/` | Banco C++ que compara un algoritmo de doble ajuste con libm, CRlibm y SLEEF; genera informe de tiempo y ULP. El propio README advierte que una auditoría formal debe incorporar MPFR. |
| Aleatorovix | `web per llançar mobles/` | Web estática y motor JavaScript sin `Math.random`; también conserva implementación Python/C de referencia. No se ha evaluado aún su calidad estadística. |
| ZypyZape / Quijote / Kilómetro | README raíz, `repo/` y gemelos | Simulaciones y propuestas de inercia sintética. Las cifras y supuestas recuperaciones de energía son afirmaciones del corpus, no validaciones independientes; requieren revisión física, balances energéticos y reproducción de simulaciones antes de tratarlas como resultados. |
| Convergencia geométrica | `encriptacionGeometrica/` y puentes AntiPC | Material educativo/experimental. No hay evidencia revisada de seguridad criptográfica; no debe presentarse como criptografía post-cuántica segura. |

El `README.md` de `repo/` mezcla documentación con una conversación/exportación extensa. Se usa para localizar conceptos y rutas, pero no como especificación canónica ni como prueba de afirmaciones científicas.

### Geometría y Aleatorovix

`encriptacionGeometrica/` es un laboratorio con muchas variantes exportadas de
chat, no una única implementación canónica. Hay puentes C de convergencia
geométrica en `repo/antipc/src/antipc/native/antipc_core/src/`, pero no se ha
encontrado una especificación de seguridad, análisis criptográfico o suite de
vectores que permita llamarlo cifrado seguro.

Aleatorovix tiene un núcleo Python y un port JavaScript. En navegador usa
`crypto.getRandomValues` cuando está disponible; en el fallback se deriva de
tiempo, jitter y estado local predecible. Sus salidas se someten además a
operaciones modulares. Por ello puede usarse como comportamiento visual o
experimental, pero no como CSPRNG ni como sustituto de `crypto.getRandomValues`
para claves, nonces, sorteos de seguridad o criptografía.

### Gemelos de energía: revisión adicional

* `VMA_mates_rescat_2026/09_gemelos_zypyzape_ultims/` contiene las variantes
  de gemelo v9–v10 y el gemelo ZypyZape/Quijote v4.8; `10_inicio...` aporta
  contexto y un ejemplo XFI. Son modelos Python/NumPy de dinámica de rotor,
  par aerodinámico, pitch y acoplamiento, no resultados de banco físico.
* `10_inicio.../02_MATH_QUIJOTE_3vs7.txt` demuestra algebraicamente que, para
  masas sinusoidales equiespaciadas, la suma de `r_k²` es constante. También
  calcula la energía cinética asociada a cambiar el momento de inercia. La
  utilidad práctica sigue dependiendo de energía/potencia del actuador,
  límites estructurales, aerodinámica y control; el propio texto distingue
  potencia transitoria de energía almacenada.
* `kilometro_sim/check_gemini_y_fisica.py` y `RESUMEN_CHEQUEO.txt` son una
  corrección interna importante: detectan que los modelos abiertos podían
  aparentar excedente al cobrar energía potencial sin rearmarla. El modelo de
  vueltas completas usa `PE = m·g·R·sin(phi)` y reporta balance mecánico;
  concluye que Kilómetro es buffer/convertidor con pérdidas, no generación
  neta de un campo gravitatorio conservativo. `README_DOBLE_KM.md` aplica la
  misma restricción al acoplamiento de dos unidades por perneo.

## Evidencia pendiente de reproducir

Los documentos registran, entre otras cifras, ~31 % de operaciones ALU evitadas
en una demo de reutilización y 62,1 % menos tiempo en un benchmark UDP frente a
una arquitectura de referencia. Son resultados históricos declarados en
`05-docs-tecnicos/RESULTADOS_BENCHMARK.txt`; faltan todavía ejecución limpia,
comando exacto, versión de código y resultados actuales para elevarlos a
evidencia reproducida.

## Próximo tramo del inventario

1. Determinar qué copia de cada proyecto es fuente activa, referencia o backup.
2. Leer y ejecutar de forma controlada las pruebas de `repo/antipc` y `vma`.
3. Mapear `VMA_mates_rescat_2026`, en especial `03_mdc`, contra las copias
   históricas L5.
4. Construir un índice de duplicados para reducir el corpus a fuentes únicas.
