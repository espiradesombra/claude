# 01_cribas

**Tres cribas VMA** + implementación C OpenMP.

## Qué hay
- `cribas.py` — Desmemoriada, Modular6k (**fase 2p/4p corregida 2026-07-23**), Híbrida
- `criba_modular.py` / `criba_hibrida.py` — scripts cortos
- `anexoF_criba6kpm1_openmp.c` — paralelo OpenMP
- Textos: desmemoriada, cotas y estructuras

## Probar
```bat
python -c "from cribas import comparar_cribas; comparar_cribas(500)"
```
