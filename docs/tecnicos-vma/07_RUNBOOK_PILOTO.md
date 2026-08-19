# Runbook de piloto  
## Checklist operativo: software → tanque Kilómetro → 2 turbinas ZypyZape

**Versión:** 1.0  
**Uso:** hoja de trabajo en laboratorio / reunión de arranque de piloto.  
**No sustituye:** licencias OEM, grid code, ni HSE de obra.  
**Documentos base:** 01 (ZZ), 03 (Kilómetro), 05 (límites), 06 (AntiPC/K3).

---

## 0. Reglas de oro (leer antes de gastar dinero)

| # | Regla |
|---|--------|
| R1 | **Nunca** publicar \(W_{net}>0\) en ciclo cerrado sin rearmar PE / stock de lastres. |
| R2 | Contar **siempre** el reset (lift lastre, actuador Quijote, motor de rearme). |
| R3 | Software y tanque se validan **antes** de pedir acceso a turbinas reales. |
| R4 | Quijote en pala **no** entra en el piloto mínimo (fase 2 opcional). |
| R5 | Uso **civil** únicamente (`02_USO_CIVIL.txt`). |
| R6 | Manifiesto SHA-256 del pack de software del piloto = línea base de auditoría. |

**Orden obligatorio de fases**

```text
F0  Software + hashes + gemelos     (días)
F1  Tanque Kilómetro MVP            (semanas)
F2  Banco rotativo / emulador ZZ    (semanas–meses)   [opcional si hay banco]
F3  Piloto 1–2 turbinas reales      (meses; con OEM/TSO)
F4  Parque / servicios ancilares    (solo tras F3 OK)
```

No saltar de F0 a F3.

---

## F0 — Software e integridad (día 0 / 15–60 min)

### F0.1 Preparación

| ☐ | Tarea | Cómo / criterio pass |
|---|--------|----------------------|
| ☐ | Python ≥ 3.10, `numpy`, `matplotlib` | `python --version` |
| ☐ | Rutas del pack en USB o `33x1\` | Índice `00_INDICE.md` presente |
| ☐ | Leer doc **05** (límites) y **06** (AntiPC) | Firma mental: sin overunity / sin RSA roto |
| ☐ | Regenerar manifiesto | `python generar_manifiesto.py` → JSON nuevo |
| ☐ | Guardar `MANIFIESTO_HASHES.json` + fecha | Baseline del piloto |

### F0.2 Demos mínimas (pass/fail)

| ☐ | Demo | Comando orientativo | Pass |
|---|------|---------------------|------|
| ☐ | Hash K3 / CLI | `python cli.py hash --text "Hola"` o equivalente AntiPC | Sale huella estable |
| ☐ | MDC toy | `mdc analyze --n 1147` o `143` | Factores correctos |
| ☐ | Red UDP 4 hubs | `network demo --hubs 4 --duration 2` | 0 drops; log JSON |
| ☐ | Benchmark A vs B o A→E | `03_benchmark_udp.bat` / `benchmark.py` | Ratio documentada en **este PC** |
| ☐ | 1ª ley Kilómetro | `python kilometro_sim\check_gemini_y_fisica.py` | `overunity_elec: false` |
| ☐ | Lastre sim | `python kilometro_sim\sim_enjambre_pesos_impar.py` | STOP al vaciar stock; saldo con lift &lt; 0 |
| ☐ | Gemelo ZZ | `python gemelo_v94.py` o twin v4.8 | PNGs generados; sin crash |
| ☐ | (opc.) XFI N=3 | `gemelo_xfi_avion_4turbinas.py --n 3` | Solo lab conceptual |

### F0.3 Entregables F0

- [ ] Carpeta `piloto_F0_YYYYMMDD\` con: manifiesto, logs JSON, capturas de bench, `RESUMEN_CHEQUEO.txt`  
- [ ] Una página: *“Qué medimos / qué no afirmamos”* (copiar tabla doc 05 §2)  
- [ ] **Go / No-Go** a F1: solo **Go** si F0.2 sin fallos críticos y manifiesto guardado  

---

## F1 — Tanque Kilómetro (MVP taller)

Referencias: `kilometro_sim\PROTOTIPO_TANQUE_PLANOS_Y_BOM.md`, `v1_1_ESP32\`, doc 03.

### F1.1 Hardware mínimo

| ☐ | Elemento | Spec MVP |
|---|----------|----------|
| ☐ | Tanque / columna | Δh útil ~0,8 m; agua 15–25 °C |
| ☐ | Módulo (objeto) | m_base 1,5–2,5 kg; guía 2 varillas |
| ☐ | Config flotación | **3** masas → flota; **+1 lastre** → se hunde |
| ☐ | Lastre transferible | 0,5–1,0 kg; 2 pernos (R y O) |
| ☐ | Make-before-break | Nunca estado 0/0 (peso suelto) |
| ☐ | Lift de recarga | **Manual** en v1 |
| ☐ | Medida | Báscula, cronómetro; opcional célula de carga / ESP32 log |
| ☐ | ESP32 v1.1 (opc.) | Solenoides O / R_alta / R_baja + FC_ALTA / FC_BAJA |

### F1.2 Calibración estática (antes de ciclar)

| ☐ | Prueba | Pass |
|---|--------|------|
| ☐ | C1: n=3 en ALTA **no se hunde** en reposo | Sí |
| ☐ | C2: n=4 **se hunde** sin empujón | Sí |
| ☐ | Margen 3→4 | +8–15 % densidad aparente documentado |
| ☐ | Seguridad: topes, no atrapamiento de manos | Checklist HSE local |

### F1.3 Ciclo de prueba (20 ciclos)

```text
1. ALTA, n=3, flota
2. Enganche lastre → n=4
3. Bajada guiada Δh
4. Entrega lastre en BAJA → n=3
5. Subida (boyancia / guía)
6. Operador devuelve lastre a ALTA (recarga batería)
7. Log: t_ciclo, caídas, fallos de perno
```

| ☐ | Criterio | Pass |
|---|----------|------|
| ☐ | C3: 0 caídas de lastre en enganche ALTA (20 ciclos) | 0 |
| ☐ | C4: 0 caídas en entrega BAJA (20 ciclos) | 0 |
| ☐ | C5: t conmutación pernos medido | tabla |
| ☐ | C6: \(E_{bajada}\) estimada vs coste de subir lastre a mano | \(E_{reset} \ge E_{util}\) (esperable) |
| ☐ | C7: **no** se reporta “motor perpetuo” | texto OK |

### F1.4 Contabilidad (plantilla)

Por ciclo y por serie de N ciclos:

| Magnitud | Símbolo | Cómo |
|----------|---------|------|
| PE lastre | \(m_L g \Delta h\) | cálculo |
| Energía recuperable estimada bajada | \(E_{down}\) | freno / estimación |
| Coste reset manual/eléctrico | \(E_{up}\) | medida o cota superior |
| Pernos | \(E_{pin}\) | nº conmutaciones × energía |
| Saldo serie con reset | \(E_{down}-E_{up}-E_{pin}\) | debe ser **≤ 0** en media |

### F1.5 Entregables F1

- [ ] Vídeo de 1 ciclo + hoja de 20 ciclos  
- [ ] BOM real vs CSV del repo  
- [ ] Log Serial ESP32 (si aplica)  
- [ ] **Go / No-Go** a F2/F3: Go solo si C1–C4 y contabilidad no miente  

---

## F2 — Banco / emulador ZypyZape (recomendado antes de turbina)

Si no hay banco, documentar el salto de riesgo y pasar a F3 solo con socio OEM.

| ☐ | Tarea | Pass |
|---|--------|------|
| ☐ | 2 ejes o 2 convertidores emulados con \(J\) conocido | dinámica ACEL/FREN visible |
| ☐ | Ciclo \(f_{ciclo}\approx0{,}4\) Hz, \(P_{ZZ,frac}\) acotada (p. ej. ≤13 %) | estable 30 min |
| ☐ | Medir \(\Delta\omega\), potencia intercambiada, temperatura | log 1 Hz |
| ☐ | Gemelo digital en paralelo (error acotado) | banda acordada (p. ej. 10–15 %) |
| ☐ | Parada de emergencia y modo “solo MPPT” | &lt; 1 s a baseline |
| ☐ | AntiPC/telemetría opcional: traza de eventos | JSON |

**No-Go F2:** oscilaciones no amortiguadas, sobrecorriente, o desacuerdo gemelo &gt; umbral sin explicación.

---

## F3 — Piloto 1–2 turbinas reales (ZypyZape software-only)

### F3.1 Pre-requisitos contractuales / técnicos

| ☐ | Ítem | Notas |
|---|------|-------|
| ☐ | Acuerdo OEM o dueño de parque (firmware / capa de control) | sin esto no hay F3 |
| ☐ | Límites de garantía AEP y fatiga por escrito | umbral p. ej. ΔAEP ≤ 0,5–1 % |
| ☐ | Acceso SCADA: \(f\), \(P\), \(\omega\), pitch, alarmas | historización ≥ 1 Hz en eventos |
| ☐ | Procedimiento de **bypass** a control nativo | 1 click / consigna |
| ☐ | Ventana de ensayo con TSO/DSO si aplica | notificar eventos de red |
| ☐ | Seguros / HSE / permisos de planta | checklist local |
| ☐ | Pack software del piloto = hash del manifiesto F0 | no “última versión del chat” |

### F3.2 Configuración mínima

| Parámetro | Valor de arranque (conservador) |
|-----------|----------------------------------|
| Turbinas | **2** adyacentes (par ACEL/FREN) |
| \(P_{ZZ,frac}\) | Empezar **5–8 %** (luego subir a ≤13 % si OK) |
| \(T_{ciclo}\) | 2,5 s (0,4 Hz) o más lento al inicio |
| Modo | **Solo evento de red** primero; ciclo continuo después |
| Quijote | **OFF** (no hardware de pala en F3 mínimo) |
| Kilómetro | **OFF** en turbina; tanque solo lab paralelo |

### F3.3 Secuencia de ensayo (días)

| Día | Acción | Criterio |
|-----|--------|----------|
| D0 | Baseline 24–72 h sin ZZ | AEP, alarmas, \(P(t)\), \(f\) local |
| D1 | ZZ en **modo sombra** (calcula pero no actúa) o par muy bajo | logs sin alarma |
| D2–D3 | ZZ activo en ventana diurna, par bajo | fatiga/temp OK |
| D4+ | Eventos de desequilibrio (naturales o acordados) | nadir/RoCoF vs baseline |
| Dn | A/B: ZZ ON vs OFF en tramos comparables de viento | informe |

### F3.4 Métricas de aceptación F3

| ☐ | Métrica | Pass (propuesta; fijar con socio) |
|---|---------|-----------------------------------|
| ☐ | Disponibilidad del control | ≥ 99 % del tiempo programado |
| ☐ | Alarmas nuevas atribuibles a ZZ | 0 críticas |
| ☐ | ΔAEP en modo continuo | dentro de contrato |
| ☐ | Mejora de nadir / RoCoF en eventos | dirección esperada **o** explicación con datos |
| ☐ | Fatiga (OEM) | dentro de curva de diseño |
| ☐ | Operador puede desactivar ZZ | verificado en ensayo |

### F3.5 Entregables F3

- [ ] Informe de 10–20 páginas: método, datos, límites (citar doc 05)  
- [ ] Series temporales y hash del firmware/capa usada  
- [ ] Recomendación **Go / No-Go** a F4 (parque)  
- [ ] Lista de bugs del gemelo vs SCADA (cerrar modelo)  

---

## F4 — Solo si F3 es Go

| ☐ | Paso |
|---|------|
| ☐ | Escalar a N=5–10 turbinas / 1 módulo |
| ☐ | Integración SCADA parque + telemetría AntiPC opcional |
| ☐ | Conversación de mercado de servicios ancilares (si existe) |
| ☐ | Quijote / Kilómetro industrial: proyectos **separados**, no mezclar con el go de ZZ software |

---

## Matriz RACI rápida (piloto)

| Actividad | VMA / inventor | OEM | Dueño parque | Lab / universidad |
|-----------|----------------|-----|--------------|-------------------|
| F0 software | R/A | C | I | C |
| F1 tanque | R/A | I | I | R (si hay taller) |
| F2 banco | R | C | I | R/A |
| F3 turbinas | C | A/R firmware | A planta | C datos |
| Informe público | R | C | A redacción | C |

R = Responsible, A = Accountable, C = Consulted, I = Informed.

---

## Plantilla Go / No-Go (copiar por fase)

```text
FASE: F_
FECHA:
RESPONSABLE:

PASS:
  [ ] Criterios técnicos C… cumplidos
  [ ] Contabilidad energética / hashes honestos
  [ ] HSE / permisos OK
  [ ] Documentación archivada en piloto_F_YYYYMMDD

RIESGOS ABIERTOS:
  1.
  2.

DECISIÓN:  GO / NO-GO / GO CONDICIONADO
CONDICIONES:
FIRMA:
```

---

## Calendario orientativo (honesto)

| Fase | Duración típica | Coste relativo |
|------|-----------------|----------------|
| F0 | 1–3 días | Muy bajo |
| F1 | 2–8 semanas | Bajo–medio (taller) |
| F2 | 1–4 meses | Medio (si hay banco) |
| F3 | 3–12 meses | Alto (acceso planta + OEM) |
| F4 | 1–3 años | Muy alto |

Acelerar F3 sin F0/F1 suele **destruir credibilidad**, no ahorrar tiempo.

---

## Checklist de reunión de arranque (1 página)

**Asistentes:** _______  
**Fecha:** _______

1. [ ] Todos leyeron **05** (límites) y aceptan no afirmar overunity.  
2. [ ] Manifiesto F0 generado y compartido.  
3. [ ] Dueño de F1 tanque y fecha de C1–C2.  
4. [ ] Nombre del contacto OEM/parque para F3 (o “aún no”).  
5. [ ] Presupuesto y qué fase se financia primero (recomendado: F0+F1).  
6. [ ] Próxima fecha de revisión Go/No-Go: _______  

---

## Referencias rápidas

| Necesitas | Documento / ruta |
|-----------|------------------|
| Física ZZ | doc 01 |
| Tanque / lastre | doc 03 + `PROTOTIPO_TANQUE_PLANOS_Y_BOM.md` |
| ESP32 | `kilometro_sim\v1_1_ESP32\` |
| Qué no vender | doc 05 |
| Bench + hashes | doc 06 |
| Comandos | `33x1\repo\33x1\03_COMANDOS.txt` |

---

*Documento 07/07 — Runbook operativo. Solo técnico. Sin promesas de mercado.*
