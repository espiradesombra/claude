# SIMULACIÓN KILÓMETRO v1 — PARÁMETROS BASE
**Víctor Manzanares Alberola — Mayo 2026**

---

## GEOMETRÍA

| Parámetro | Valor | Justificación |
|-----------|-------|--------------|
| **Forma recorrido** | Carril circular (tiovivo rotatorio) | Más simple que tubo espiral, conserva la idea de "girar" |
| **Radio de rotación** | R = 5 m | Escala industrial pequeña (no es ni prototipo ni megamáquina) |
| **Longitud del carril** | L = 2 m | Espacio para que la masa se deslice sin "caer" |
| **Número de masas** | N = 3 | Coincide con lógica 3 pales de Quijote (desfase 120°) |
| **Masa por bola** | m = 500 kg | Densa (acero) pero manejable mecánicamente |

---

## CICLO DE OPERACIÓN

| Parámetro | Valor | Justificación |
|-----------|-------|--------------|
| **Recorrido radial mín/máx** | r_min = 3 m, r_max = 7 m | 4 metros de recorrido (razonable para 500kg) |
| **Velocidad angular** | ω_rot = 1.0 rad/s | 0.159 rpm ≈ 1 revolucción cada 6.28 segundos (lento, controlable) |
| **Período de ciclo** | T_ciclo = 2π/ω = 6.28 s | 1 ciclo cada ~6 segundos |
| **Ciclos/hora** | 570 ciclos/hora | Operación continua |

---

## FUERZAS Y CONTROL

| Parámetro | Valor | Justificación |
|-----------|-------|--------------|
| **Aceleración gravitatoria** | g = 9.81 m/s² | Estándar |
| **Densidad fluido (agua)** | ρ = 1000 kg/m³ | Simulamos "bajo el mar" |
| **Coef. fricción carril/masa** | μ_k = 0.05 | Cojinetes de rodillo (bajo) |
| **Factor de amortiguamiento** | c = 100 N·s/m | Viscosidad del fluido + fricción |
| **Control de velocidad** | v_slide(t) = A·sin(ω·t + φ_i) | Sinusoidal coordinada por masa i |
| **Amplitud de deslizamiento** | A_slide = 0.5 m/s | Máximo 0.5 m/s en carril |

---

## ENERGÍA Y POTENCIA

| Parámetro | Valor | Justificación |
|-----------|-------|--------------|
| **Potencia nominal esperada** | P_nom = 50 kW | Módulo pequeño (comparar con Quijote 2.5 MW × 10 turbinas = 25 MW) |
| **Eficiencia simulada** | η_sim_target = 60-70% | Realista con fricción + amortiguamiento |
| **Energía almacenada / ciclo** | E_ciclo = ? | **A calcular** |

---

## ENSAMBLE ZYPYZAPE

| Parámetro | Valor | Justificación |
|-----------|-------|--------------|
| **Número de Kilómetros** | N_KM = 5 | Módulo pequeño (como 5 aerogeneradores) |
| **Intercambio de energía** | Ciclo: 1 genera, otros 4 consumen/cargan | ZypyZape → repartición temporal |
| **Frecuencia de intercambio** | f_intercambio = 1/6.28 Hz ≈ 0.159 Hz | Sincronizado con ciclo Kilómetro |

---

## SIMULACIÓN

| Parámetro | Valor | Justificación |
|-----------|-------|--------------|
| **Tiempo total simulado** | T_sim = 1 día (86400 s) | ≈ 13,760 ciclos |
| **Paso de integración** | dt = 0.01 s | Exactitud suficiente para dinámicas de ~1 Hz |
| **Número de pasos** | N_steps = 8,640,000 | 1 día @ 0.01 s/paso |

---

## SALIDAS ESPERADAS

1. **Gráfica energía**: E_entrada(t), E_salida(t), E_red(t), eficiencia_neta(t)
2. **Gráfica posición**: r_i(t) para cada masa i
3. **Gráfica velocidad**: ω_rot(t), v_slide(t)
4. **Gráfica fuerza**: F_centrífuga(t), F_fricción(t), F_neta(t)
5. **Gráfica balance**: P_generat(t) vs P_consumida(t)
6. **Estadísticas**: η_media, η_min, η_max, E_total_generada, E_total_consumida

---

**PRÓXIMO PASO:** Código MATLAB/Octave que implemente esto.
