# docs — física VMA (ZypyZape · Quijote · Kilómetro)

## Lectura principal

| Archivo | Qué es |
|---------|--------|
| [EXPLICACION_FISICA_ZYPYZAPE_QUIJOTE_KILOMETRO.md](./EXPLICACION_FISICA_ZYPYZAPE_QUIJOTE_KILOMETRO.md) | Texto unificado (+ 1477) |
| [diagrama_fisica_ZZ_Quijote_Kilometro.png](./diagrama_fisica_ZZ_Quijote_Kilometro.png) | **Diagrama 1 página** |
| [gemelo_kilometro_minimo.py](./gemelo_kilometro_minimo.py) | **Gemelo Kilómetro** balance 1ª ley |
| [kilometro_minimo_balance.png](./kilometro_minimo_balance.png) | Gráfica del gemelo |

## Generar / ejecutar

```bat
cd docs
python generate_diagrama_fisica_1pagina.py
python gemelo_kilometro_minimo.py --mode drain --T 25
python gemelo_kilometro_minimo.py --mode maintain --T 25 --omega0 2
```

### Gemelo Kilómetro — modos

| mode | Comportamiento |
|------|----------------|
| `drain` | Arranca con ω alta y **vacía** la batería cinética → E_útil ≈ η_reg · dissipación |
| `maintain` | Actuador mantiene ω → hace falta **E_ext** continua |

**Resultado típico (drain):** η_paid = E_util/(E_ext+drain) ≈ η_reg ≤ 1.  
No hay motor perpetuo: se gasta Emec inicial.

## Enlaces

- [XFI](../XFI.md)
- [33×1 qué es](../33x1/00_QUE_ES_33x1.md)
- Pack rescat: [VMA_mates_rescat_2026/](../VMA_mates_rescat_2026/)
