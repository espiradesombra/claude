# Resultados toy — Quijote modes + ZZ×3 calibrado  
**2026-07-23** · escala didáctica (no NREL full)

## Quijote 3 vs 7 — coste de actuador

```bat
python gemelo_quijote_3vs7_balance.py --T 20 --mode static|phase|ball
```

| Mode | N | P_gen μ | P_act μ | **P_net (gen−act)** |
|------|---|---------|---------|---------------------|
| **static** | 3 | 308 W | **0** | **+308 W** |
| **static** | 7 | 302 W | **0** | **+302 W** |
| **phase** | 3 | 309 W | 575 W | **−265 W** |
| **phase** | 7 | 300 W | 1314 W | **−1014 W** |
| **ball** | 3 | 308 W | 1214 W | **−906 W** |
| **ball** | 7 | 302 W | 2760 W | **−2459 W** |

### Lectura

1. **static** = turbina sin mover masas → referencia limpia (P_act = 0).  
2. **phase** (maniobras discretas + histéresis) cuesta **menos** que **ball** continuo.  
3. **7 palas** en este toy **no gana** en P_net: hay **más masas** que mover → más P_act.  
4. El valor de 7 palas en el paper es **continuidad / rizado**, no “más watios gratis”.  
5. Si el control no aporta más captación aero o servicios de red valorados, el actuador **come el margen**.

PNG: `quijote_3vs7_balance.png` (último run = mode phase).

---

## Grupo ZZ×3 — H calibrada

```bat
python gemelo_grupo_zz3.py --T 35 --compare-H
```

| Métrica | ZZ ON | ZZ OFF |
|---------|-------|--------|
| f min / max | 49.86 / 50.00 Hz | (similar) |
| **nadir** post escalón +15% | **49.859 Hz** | **49.858 Hz** |
| Δnadir (ON−OFF) | **+0.001 Hz** | — |
| RoCoF max | ~0.18 Hz/s | — |
| H_ZZ mean | ~0.39 s | — |

### Lectura

1. Nadir en **décimas de Hz**, no caídas irreales a 45 Hz (calibración OK).  
2. Beneficio ZZ en este toy es **pequeño pero positivo** (+1 mHz) — coherente con twins v4.8 “casi invisible en 2 GW / marginal por módulo”.  
3. El valor de ZZ es **escalar módulos** y servicios ancilares, no un milagro de frecuencia en un solo grupo de 3.

PNG: `grupo_zz3_balance.png`

---

## Mensaje unificado

| Sistema | Qué queda en pie |
|---------|------------------|
| Kilómetro | Batería / convertidor; η_paid ≤ η_reg |
| Quijote | Útil si el **control** vale más que **P_act**; phase > ball; static = base |
| ZZ×3 | Inercia sintética y roles; mejora de nadir **marginal** por grupo pequeño |

Todo es **toy**. Piloto real = siguiente nivel de credibilidad.
