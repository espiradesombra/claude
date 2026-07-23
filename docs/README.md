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
python gemelo_quijote_3vs7_balance.py --T 25 --mode ball
python gemelo_grupo_zz3.py --T 30 --K 80
```

### Gemelos

| Script | Qué demuestra |
|--------|----------------|
| `gemelo_kilometro_minimo.py` | 1ª ley; η_paid ≤ η_reg (batería cinética) |
| `gemelo_quijote_3vs7_balance.py` | 3 vs 7: `static` / `phase` / `ball` + coste actuador |
| `gemelo_grupo_zz3.py` | ZZ×3 calibrado; `--compare-H` ON vs OFF |
| `RESULTADOS_TOY_2026-07-23.md` | **Tabla de números** de la última corrida |

| mode (Kilómetro) | Comportamiento |
|------|----------------|
| `drain` | Vacía batería cinética → E_útil ≈ η_reg · disipación |
| `maintain` | Actuador mantiene ω → E_ext continua |

**Mensaje común:** no hay motor perpetuo; actuador y fricción cuentan.

## Enlaces

- [XFI](../XFI.md)
- [33×1 qué es](../33x1/00_QUE_ES_33x1.md)
- Pack rescat: [VMA_mates_rescat_2026/](../VMA_mates_rescat_2026/)

## Meta

- [DE_ON_VE_TOT_AIXO.md](./DE_ON_VE_TOT_AIXO.md) — d'on surt tot, què deien les IAs, criteris d'honestedat

