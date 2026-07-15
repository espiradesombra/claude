# Tabla de contenido — Ecosistema VMA / 33×1

**Autor:** Víctor Manzanares Alberola (EPSA, UPV Alcoi)  
**Paquete:** `just run` — índice maestro por tipo de activo  
**Versión:** 2026-07-11

---

## Leyenda de estado

| Símbolo | Significado |
|---------|-------------|
| ✓ | Implementado / reproducible en `just run` |
| ◐ | Documentado, simulación o heurística; falta cierre formal/hardware |
| ○ | Conjetura o borrador; en investigación |
| $ | Mercado identificado (ver columna mercados) |

---

## 1. CONJETURAS

| ID | Nombre | Enunciado breve | Estado | Dónde |
|----|--------|-----------------|--------|-------|
| C01 | Intervalo mínimo dos primos | En \(I(n)=[n-\lfloor\sqrt{n+3}\rfloor-3,\,n+3]\) hay ≥2 primos si \(L(n)-m(n)\geq2\) | ◐ | `archivos-vma/cribas_cotas_vma.txt` §2–3 |
| C02 | Andrica vía salto | \(\text{SaltoMáx}(n)\approx\sqrt{n}\); coherencia asintótica con \(\sqrt{n}-\sqrt{n-\sqrt{n}}<1\) | ◐ | `cribas_cotas_vma.txt` §4 |
| C03 | Goldbach fuerte (bits) | Secuencias NOT/suma exploran patrones; relación informal con pares primos | ○ | `archivos-vma/conjeturas/`, `repe/*.xlsx` |
| C04 | Intervalo factorial \(n!\pm n\) | Entre \(n!+n\) y \(n!+\text{SaltoMáx}(n!)+n\) hay ≥2 primos \(6k\pm1\) consecutivos | ○ | `conjetura_resumen.txt` |
| C05 | Siguiente primo dado \(p\) | Heurística de predicción del primo siguiente | ○ | `cribas_cotas_vma.txt` §13, `siguiente_primo.py` |
| C06 | Densidad «cine» | Intuición: de \(n\) a \(2n\) vs de \(n\) a \(2\) — asimetría de densidad de primos | ○ | `matematica_extra/densidad_x_cine.txt` |
| C07 | Antiprimos \(A^*\) | \(AM(p,v)=(2^p+1)/3^v\) fraccionario → \(A^*=\lceil AM\rceil_3\cdot AM\) revela factorización | ○ | `archivos-vma/antiprimos.txt` |
| C08 | Hurto gravitatorio neto | Ciclo 540°+180° / 3 vs 7 pales → \(P_{net}>0\) en condiciones PAPER | ◐ | `hurto-gravitatorio/quijote_results.csv` |
| C09 | Captura 150% PIB | Ruptura de simetría temporal en campo conservativo → reordenación económica global | ○ | `33x1/00_PROMPT_DEFINITIVO_33x1.md` |
| C10 | 33 años paz racional | Adopción Quijote hace guerra 3× menos rentable que paz para rivales | ○ | `33x1/00_PROMPT_DEFINITIVO_33x1.md` §V |

---

## 2. NÚMEROS (estructuras, familias, constantes)

| ID | Objeto | Definición / rol | Dónde |
|----|--------|------------------|-------|
| N01 | Candidatos \(6k\pm1\) | Rejilla modular de primos \(>3\) | `cribas`, `criba_hibrida.py`, Aleatorovix |
| N02 | \(L(n)\) | \(\lfloor\sqrt{n+3}\rfloor+7\) — longitud de intervalo | `anexoE_L_m_script.py` |
| N03 | \(m(n)\) | \(\sum_{i=2}^{K}(\sqrt{n+3}/i!)\approx(e-2)\sqrt{n}\) — marcado compuestos | `anexoF_mrauv_calibrador.py` |
| N04 | \(K(n)\) | \(\lfloor\lfloor\sqrt{n}\rfloor\cdot 9/24\rfloor\) — tope sumandos | `anexoE_L_m_script.py` |
| N05 | Factor\_N | \(2\sin(\pi/N)\) — asimetría 3 vs 7 pales | `zypyzape-contexto/02_MATH_QUIJOTE_3vs7.txt` |
| N06 | Familias Sofí | \(L_1,L_2,L_3,L_4, L_{SG}\) — Sophie Germain modular | `estructura_sofi.py`, `cribas_cotas_vma.txt` §5 |
| N07 | Férmines modulares | Alineación \(6k\pm1\) en números de Fermat | `fermat_modular.py` |
| N08 | Primos gemelos / gaussianos | Proporción factorial entre pares | `jj.txt`, `numeros_i_numeritos_jj.txt` |
| N09 | 1477 / bandera fase | Codificación booleana qubit \((v,f)\); distancia Hamming | `33x1/1477.md` |
| N10 | 33×1 | 33 años paz ↔ 1 adhesión territorial irrevocable | `33x1/` |
| N11 | MDC restos | \(n \bmod i\) + hipotenusa — factorización geométrica | `vma-run/lib/mdc.py`, `restos_graficos.py` |
| N12 | \(\pi_6(n)\) Criva | Densidad racional iterativa \(D_n=(a\cdot p+b)/(b\cdot p)\) | `criva_teoria.txt`, `teoremas/criva.py` |

---

## 3. TEOREMAS (demostrados o formalizados)

| ID | Nombre | Resultado | Estado | Dónde |
|----|--------|----------|--------|-------|
| T01 | Amplificador fase K=3 XOR | Realimentación booleana Toffoli vs AND clásico → bandera de coherencia | ✓ | `teoremas/THEOREM_PhaseAmplifier_K3_XOR.md` |
| T02 | \(L_4\setminus L_2\) infinito | Infinitud vía TCR + Dirichlet en clases Sofí | ◐ | `cribas_cotas_vma.txt` §5 |
| T03 | Criterio dos primos | Si \(L(n)-m(n)\geq2\) ⇒ ≥2 primos en \(I(n)\) (operativo) | ◐ | `cribas_cotas_vma.txt` §2 |
| T04 | MDC factorización | \(10403=101\times103\) — método cinemático reproducible | ✓ | `RUN.bat mdc 10403` |
| T05 | AntiPC H1–H3 | Desacoplamiento + knowledge \(K(N)\) reduce carga ALU medible | ✓ | `antipc/docs/explicacion_cientifica.txt` |
| T06 | Conservación ForPkm | Ciclo cerrado gravedad: sin aporte externo \(W_{neto}\not>0\) (ideal) | ◐ | `energia/giroscopio_forpkm.txt` |
| T07 | K3 auditoría | Hash/fingerprint industrial — no criptográfico, determinista | ✓ | `vma-k3/` |

---

## 4. LEMAS (auxiliares, axiomas, principios)

| ID | Nombre | Enunciado | Ámbito | Dónde |
|----|--------|-----------|--------|-------|
| L01 | Candidato 6k±1 | Todo primo \(p>3\) ≡ \(1\) o \(5\) (mod 6) | Número | `cribas_cotas_vma.txt` §1 |
| L02 | Salto máximo | \(\text{SaltoMáx}(n)\lesssim\sqrt{n}\) (cota superior) | Número | `cribas_cotas_vma.txt` §2 |
| L03 | Contracción Criva | Error reduce ~50% por iteración con parámetro \(s\) | Número | `criva_teoria.txt` |
| L04 | Ley consulta AntiPC | Operar solo si knowledge no resuelve | Protocolo | `gptcomputing/CONCEPTO.txt` |
| L05 | Ley inmutabilidad | Reference `frozen=True` — conocimiento publicado no muta | Protocolo | `gptcomputing/MAPA_*.txt` |
| L06 | Ley activación | EventBus por eventos, sin polling | Protocolo | `antipc/runtime/event_bus.py` |
| L07 | Axioma 11 ciclo vida | Runtime por estado (Epoch, Tick), no reloj físico | Protocolo | `gptcomputing.txt` |
| L08 | Axioma 17 consenso | Conocimiento consolidado = evaluadores independientes coinciden | Protocolo | `gptcomputing.txt` |
| L09 | Factor fase Quijote | 1 pala expande / \((N-1)\) contraen → asimetría temporal | Renovable | `33x1/00_PROMPT_DEFINITIVO_33x1.md` |
| L10 | Entropía Aleatorovix | Bits físicos (reloj, pila, heap) alimentan diccionario mutante | Método | `aleatorovix/CONCEPTO.txt` |
| L11 | Knowledge Value | \((Reuse\times Cost\times Trust)/Storage\) | Protocolo | `gptcomputing.txt` v0.0.072 |
| L12 | El uno molesta | Paquete comercial indivisible 33×1 ≠ licencia parcial | Mercado | `33x1/PAPER_33x1_IRONIA.txt` |

---

## 5. MÉTODOS

| ID | Nombre | Qué hace | Implementación | Dónde |
|----|--------|----------|----------------|-------|
| M01 | MDC | Factorización por saltos hipotenusa / cinemática | ✓ Python | `vma-run/lib/mdc.py` |
| M02 | Criva racional | Densidad primos por iteración primorial top-down | ✓ Python | `teoremas/criva.py`, `vma-run/lib/criva.py` |
| M03 | Criba híbrida 6k±1 | Marcado + salto \(2p/4p\) alternado | ✓ Python | `archivos-vma/codigo/criba_hibrida.py` |
| M04 | Criba desmemoriada | Patrón binario distribuido; `memset` suicida | ✓ Py + C | `aleatorovix/` |
| M05 | Criba OpenMP | 6k±1 paralela, triángulo superior | ✓ C | `anexoF_criba6kpm1_openmp.c` |
| M06 | MRAUV | Densidad por tramos; calibrador \(V_0,a_0\) | ✓ Python | `mrauv.py`, `anexoF_mrauv_calibrador.py` |
| M07 | Restos + hipotenusa | Visual / factorización modular gráfica | ✓ Python | `restos_graficos.py`, `gen_fotos.py` |
| M08 | K3 industrial | Modos usual/propio-b; auditoría memoria | ✓ pip | `vma-k3/` |
| M09 | AntiPC runtime | Flow kernel + plugins + UDP fabric | ✓ Python | `antipc/` |
| M10 | Gemelo Quijote | NREL 5MW + buffer hidráulico + \(r_q\) | ✓ Python | `gemelos/gemelo_v6_fisic.py` |
| M11 | Red ZypyZape | Swing equation + inercia sintética v4.8 | ✓ Python | `zypyzape_twin_v4_8_quijote.py` |
| M12 | Hurto paper rules | Tabla NREL, modos PAPER/SINUS/CONTROL | ✓ Python | `hurto-gravitatorio/` |
| M13 | ForPkm | Kilómetro asíncrono 2–3 módulos desfasados | ◐ texto | `archivos-vma/energia/` |
| M14 | Phase amplifier K=3 | XOR feedback 6 ciclos, flag Hamming | ✓ | `vma-run/lib/k3_phase.py` |
| M15 | Aleatorovix Lila | Criba / salto 6k±1 / resonancia ZypyZape | ◐ C+Py | `aleatorovix/organismo_lila_v99.c` |
| M16 | Newton rápido / oráculo | Búsqueda con oráculo (Libro6) | ○ | `claude-main-extract` (no empaquetado) |
| M17 | LoveEarthHacker LaTeX | Libro métodos simples unificados | ◐ TeX | `libro-metodos-simples/` |

---

## 6. PROTOCOLOS

| ID | Nombre | Capas | Estado | Dónde |
|----|--------|-------|--------|-------|
| P01 | AntiPC v0.1.0 | Kernel → Knowledge → Plugin (UPS) → EventBus → UDP | ✓ | `antipc/`, `gptcomputing/MAPA_*.txt` |
| P02 | Arquitecturas A–E | Convencional / lock-free / hubs / Grafcet / UDP real | ✓ bench | `antipc/architectures.py` |
| P03 | K3 plugin chain | HASH → dedup → file → MDC plugins | ✓ | `antipc/plugins/` |
| P04 | Phase-Boolean 1477 | \((v,f)\) + shift Gray + Toffoli + feedback 6 ciclos | ✓ | `33x1/1477.md`, teorema T01 |
| P05 | UDP fabric AntiPC | Hub 1..N → maestro; benchmark ALU | ✓ | `udp_benchmark.py` |
| P06 | Adhesión 33×1 | Transferencia PI territorial perpetua ↔ 33 años paz | ○ borrador | `33x1/`, pendiente contrato |
| P07 | VMA-RUN demo | Un comando → 6 piezas ecosistema | ✓ | `RUN.bat`, `vma_run.py` |
| P08 | Criba desmemoriada red | Workers con memoria propia; mensajes entre nodos | ◐ | `criba_desmemoriada_guion.txt` |
| P09 | Elewit dossier | Memo 1 pág + formulario + guion 3 min | ✓ | `zypyzape-contexto/dossier-elewit/` |
| P10 | Benchmark o research | Nada entra al núcleo sin medición | ✓ norma | `gptcomputing.txt` final |

---

## 7. RENOVABLE (energía, gemelos, inercia)

| ID | Activo | Parámetros / resultado | Estado | Dónde |
|----|--------|------------------------|--------|-------|
| R01 | Turbina NREL 5MW | Base física gemelos | ✓ | `gemelos/`, contexto |
| R02 | Quijote 3 pales PAPER | \(P_{net}=+9\) kW, \(\eta_{hurto}=3.0\times\) | ◐ sim | `quijote_results.csv` |
| R03 | Quijote 7 pales PAPER | \(P_{net}=+31.1\) kW, \(\eta_{hurto}=3.7\times\) | ◐ sim | `quijote_results.csv` |
| R04 | Modo CONTROL | Centrífuga práctica — aún negativo | ○ brecha | `hurto-gravitatorio/` |
| R05 | Buffer hidráulico | \(r_q\), válvulas, estados v6/v942 | ✓ sim | `gemelo_v6_fisic.py` |
| R06 | ZypyZape inercia | Factor 3 vs 7, ball vs estático | ✓ | `02_MATH_QUIJOTE_3vs7.txt` |
| R07 | Kuramoto / red | Sincronización parque | ✓ sim | `zypyzape_twin_v4_8_quijote.py` |
| R08 | Helicoidal 540°+180° | Batería virtual + modo ZypyZape | ◐ doc | `informe_helicoidal_zypyzape.txt` |
| R09 | ForPkm → batería | Misma máquina como almacén cinético | ◐ | `energia/pacopi_forpkm.txt` |
| R10 | Servicios red | FCR/FFR/inercia ~500–900 €/MW·año EU | $ ref | `resumen de dinero...txt` |
| R11 | Pitch + hidráulica | Control actuador auditado | ◐ | `quijote_actuador.py.txt` |

---

## 8. MERCADOS_AFECTADOS

| Mercado | Activos VMA que encajan | TRL hoy | Ticket orientativo | Bloqueo principal |
|---------|-------------------------|---------|-------------------|-------------------|
| **Red eléctrica / TSO** | R02–R07, P09, ZypyZape | 4 sim | Piloto 5–20M € | Certificación IEC, CONTROL negativo |
| **OEM aerogenerador** | Gemelos, hurto, 3vs7 | 4 sim | Licencia / joint | Validación hardware |
| **Servicios ancilares** | Inercia sintética, FCR | 3–4 | 500–900 €/MW·año | Demostrar en parque real |
| **Industrial IoT / auditoría** | K3, M08, T07 | 4 | 20k–150k €/año | 1 cliente piloto |
| **Edge / DevOps** | AntiPC, P01–P05 | 4–5 | Niche B2B | Estandarizar producto |
| **Ciber / fingerprint** | K3, phase 1477 | 3 | Research | No vender como post-RSA aún |
| **Math software / HPC** | M01–M07, Criva | 3 | Académico / nicho | Peer review |
| **Deeptech opción** | Portfolio completo | — | 500k–3M € opción | Due diligence |
| **Soberanía territorial** | P06, 33×1, 1477 | 0 contrato | Paquete político | Redactar adhesión |
| **REE / Elewit** | Dossier, CSV hurto | 4 sim | Credibilidad | Expectativas perpetuum |
| **UPV / nube didáctica** | Separado de PI | — | Social | No confundir con producto |

---

## 9. MAPA RÁPIDO — qué leer primero

| Si te interesa… | Empieza por… |
|----------------|--------------|
| Números / conjeturas | `archivos-vma/cribas_cotas_vma.txt` |
| Teorema verificado | `teoremas/THEOREM_PhaseAmplifier_K3_XOR.md` |
| Demo 10 s | `RUN.bat` |
| Renovable / REE | `zypyzape-contexto/01_CONTEXT_ZYPYZAPE.txt` + CSV hurto |
| Protocolo software | `gptcomputing/CONCEPTO.txt` → `antipc/INICIO.bat` |
| Mercado / dinero | `zypyzape-contexto/resumen de dinero...txt` |
| 33×1 / territorial | `33x1/00_PROMPT_DEFINITIVO_33x1.md` |
| Subida GitHub/Elewit | `README.md` |

---

## 10. Pendientes (huecos del índice)

- [ ] Contrato adhesión territorial 33×1 (P06)
- [ ] Empaquetar Goldbach (55 archivos, C03)
- [ ] Cerrar R04 modo CONTROL hurto
- [ ] Peer review Criva / conjeturas C01–C07
- [ ] `kilometre;batería` (38 archivos) → `archivos-vma/energia/`
- [ ] Golden tests K3 Python↔C

---

*Índice vivo — actualizar al añadir carpetas a `just run`. Ver `MANIFEST.txt`.*