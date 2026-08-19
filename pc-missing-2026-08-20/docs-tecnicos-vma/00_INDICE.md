# Dossier técnico VMA — Paquete de presentación (v1.3)

**Alcance:** documentación técnica de bloques de energía, control e informática civil.  
**Enfoque:** solo ingeniería, física, matemática y estado de madurez.  
**Fecha de compilación:** 2026-08-11 (v1.3: + doc 07 runbook piloto)  
**Autor conceptual:** Víctor Manzanares Alberola (EPSA / UPV Alcoi)  
**Corpus de referencia en disco:** `Desktop\33x1\`, `Desktop\kilometro_sim\`, monorepo enlazado.

---

## Documentos incluidos

| # | Archivo | Bloque | Madurez (honesta) |
|---|---------|--------|-------------------|
| 01 | `01_ZYPYZAPE_sistema_control_inercia.md` | Control de parques / batería cinética virtual | TRL ~4–5 simulación NREL; sin piloto de campo archivado |
| 02 | `02_QUIJOTE_pesos_moviles_aspas.md` | Masas móviles en palas; \(J(t)\) | Formalismo + gemelos; hurto gravitatorio **neto = 0** en ciclo cerrado |
| 03 | `03_KILOMETRE_gravedad_flotabilidad.md` | Kilómetro / lastre / perneo / flotación | Sims 1ª ley + diseño tanque/ESP32; **no** overunity |
| 04 | `04_AVANCES_MATEMATICOS_MDC_geometria_Pi_E.md` | MDC, geometría, π/e, cifrado de fase | Código ejecutable; toy en factorización; uso civil |
| 05 | `05_LOGROS_LIMITES_TRANSICION_Y_DESNUCLEARIZACION.md` | Qué se puede / no se puede | Tabla demostrado/heurístico/refutado |
| 06 | `06_ANTIPC_K3_INTEGRIDAD_Y_ARQUITECTURA.md` | AntiPC, K3, manifiesto, L0–L4, XFI breve | TRL software más alto; demoable hoy |
| 07 | `07_RUNBOOK_PILOTO.md` | **Checklist F0→F4:** software, tanque, 2 turbinas | Operativo; no sustituye OEM/HSE |

---

## Lectura recomendada según perfil

| Perfil | Orden |
|--------|-------|
| Empezar a ejecutar | **07** → 06 → 05 → 03 → 01 |
| Revisor técnico exigente | **05** → **06** → 03 → 01 → 02 → 04 → 07 |
| Integrador software | **06** → 04 → 07 (F0) |
| Ingeniería eólica / red | 01 → 02 → 05 → 07 (F3) → 03 |
| Prototipo físico / tanque | 03 → **07** (F1) → 05 |
| Matemática | 04 → 06 → 05 |
| Transición / nuclear apagada | **05** → 01 → 03 → 06 |

---

## Principios de honestidad del dossier (v1.3)

1. Gravedad/flotación ideales: **ciclo cerrado → \(W=0\)**.  
2. Kilómetro = batería/buffer; **overunity = false**.  
3. ZypyZape = control de inercia; efecto multi-GW pequeño por módulo.  
4. Quijote = \(J(t)\) local; no gana RoCoF ms vs BESS.  
5. Desnuclearización = MWh + servicios; packs ayudan a **servicios**, no baseload mágico.  
6. AntiPC/K3 = medible y hasheable; no criptografía de Estado.  
7. Piloto: **F0 → F1 → (F2) → F3**; no saltar a turbinas sin software ni tanque (doc 07).  
8. Uso civil: `02_USO_CIVIL.txt`. Cruzar números con **05** + `LIMITES_HONESTOS.txt`.

---

## Fuentes clave

| Fuente | Aporte |
|--------|--------|
| `Desktop\kilometro_sim\` | 1ª ley, lastre 3/4, tanque BOM, ESP32 |
| `33x1\repo\33x1\` | Pack auditable, uso civil, comandos |
| `33x1\repo\teoremasgrok\` | TRL por idea, límites honestos |
| `33x1\antipc\` + `05-docs-tecnicos\` | Bench AntiPC, arquitectura |

---

## Cómo reproducir chequeos clave

```bat
cd C:\Users\cuent\Desktop\kilometro_sim
python check_gemini_y_fisica.py
python sim_enjambre_pesos_impar.py

cd C:\Users\cuent\Desktop\33x1\repo\33x1
python generar_manifiesto.py
```

Detalle de fases y checklists: **`07_RUNBOOK_PILOTO.md`**.

---

*Dossier v1.3 — solo contenido técnico.*
