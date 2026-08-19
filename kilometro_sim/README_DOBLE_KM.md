# Doble Kilómetro 90° + Perneo

## Regla de regeneración (pedido explícito)

| Fase | Condición | Actuación |
|------|-----------|-----------|
| **Desfavorable** | `cos(φ) ≥ 0` → gravedad se opone a `+ω` | **Motor o libre. REGEN OFF** |
| **Favorable** | `cos(φ) < 0` → gravedad asiste | Regeneración ON (+ seguimiento de carga) |

Código: `control_module()` en `doble_km_90_perneo.py`.  
Verificación automática: 0 regeneraciones en muestras desfavorables.

## Qué optimiza (y qué no)

- **No** busca `W_net > 0` (1.ª ley: imposible en ciclo cerrado con campo conservativo).
- **Sí** busca:
  - suavizar entrega eléctrica / reducir error ante picos,
  - usar **perneo/desperneo** como supercondensador mecánico,
  - métricas `P_pico/P_media`, error en pulso, `ω_min` durante descarga.

## Perfil de carga por defecto

- Base: 80 W  
- Pulso descarga: **2.5 kW × 0.8 s = 2000 J**  
- Pulso absorción: **−1.8 kW × 0.6 s**

Cambia `LoadPulse` en el script para otros kW/duración.

## Cómo ejecutar

```bash
python kilometro_sim/doble_km_90_perneo.py
```

Salidas:

- `doble_km_90_perneo.png`
- `doble_km_90_perneo_resultados.json`

## Nota sobre 90°

Con dos masas a 90°:

```
τ_g,total = −m g R √2 cos(φ + π/4)
```

El peak-to-peak de par gravitatorio **crece ×√2** respecto al mono.  
La utilidad del 90° es **complementariedad de fases** (uno en favorable mientras el otro está en desfavorable a menudo) y el **embrague de perneo**, no la cancelación del armónico (eso es **3 módulos a 120°**).

## Honestidad de escala

Si te dijeron que Quijote es viable con un molino gigante pero no con molinos/actuadores “de la Tierra”, el cuello real es:

`par_perturbación (viento/gravedad)  vs  par_actuador disponible`

Kilómetro **no rompe** ese cuello con sobreunidad; aporta **buffer de potencia e inercia**. Eso sí es ingeniería útil.
