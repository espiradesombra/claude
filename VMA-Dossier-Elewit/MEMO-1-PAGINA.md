# ZypyZape — Memoria ejecutiva (1 página)

**Solicitante:** Víctor Manzanares Alberola · NIF 48598033E  
**Contacto:** victormanzanaresalberola@gmail.com · +34 674 941 533  
**Área:** EPSA / investigación independiente · El Palomar (Alicante)

---

## Qué es

**ZypyZape** es un algoritmo de control para **parques eólicos existentes** que coordina el par de las turbinas (roles CAPT / ACEL / FREN) y aumenta la **inercia sintética** de la red sin baterías de litio ni hardware nuevo en la **versión base** (solo firmware / capa de control).

- Cumplimiento orientado a **ENTSO-E NC RfG** (inercia sintética).
- Gemelo digital calibrado frente referencia **NREL 5 MW** (dominio público).
- Topología modular: pares de turbinas + acoplamiento tipo **Kuramoto**.

## Qué demostramos hoy (reproducible en este ZIP)

| Entregable | Archivo |
|------------|---------|
| Gemelo completo (paper) | `gemelo_v94.py` |
| Gemelo validación NREL | `zypyzape_twin_v4_8_quijote.py` |
| Contexto y parámetros | `01_CONTEXT_ZYPYZAPE.txt` |
| Gráficos | `zypyzape_v4_validacion.png`, `zypyzape_v4_dinamica.png`, `gemelo_v8_final.png` |
| Propuesta ampliada | `ZypyZape_Propuesta_Elewit_REE_2026.docx` |

**Ejecución:** `pip install numpy matplotlib` → `python gemelo_v94.py` (genera figuras en carpeta).

## Resultados simulación (honestos)

- Red modelo 2 GW, perturbación −100 MW: mejora de **nadir** de frecuencia ~**+0,002 Hz por módulo** (10 turbinas) — efecto **lineal** con número de módulos.
- **Quijote** (masas en palas): valor principal en **dinámica local** y microredes; en red 2 GW el aporte es **marginal** (documentado en `01_CONTEXT`).
- Ganancia neta de captación eólica: **no demostrada en campo** — requiere piloto.

## Qué NO es

- No sustituye un Megapack electroquímico 1:1.
- No rompe RSA ni incluye “paquete geopolítico 33×1” en esta solicitud.
- No incluye certificación IEC ni acceso OEM — eso requiere **socio industrial**.

## Qué pedimos a Elewit / Redeia

1. **Reunión de valoración técnica** (60 min) con ingeniero de red / eólica.
2. **Carta de encaje** con programa *venture client* o línea de innovación REE, si procede.
3. **Piloto en 1–2 turbinas** (financiación y acceso PLC por el socio — no por el inventor).

**No pedimos inversión en efectivo al inventor.**

## Modelo a largo plazo (solo tras piloto)

- Estudio / consultoría inicial.
- Licencia por MW instalado o revenue-share en servicios ancilares (inercia sintética).
- Fase 2 opcional: integración **Quijote** (hardware en palas) — fuera de alcance de esta fase.

## Inventor

Trabajo de años en simulación y control; afiliación académica EPSA-UPV Alcoi. Código y memoria adjuntos para revisión por terceros.

— Víctor Manzanares Alberola · 2026