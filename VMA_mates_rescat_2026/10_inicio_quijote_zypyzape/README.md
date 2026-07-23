# 10_inicio_quijote_zypyzape

**Inicios conceptuales** Quijote / ZypyZape + **gemelo XFI** (avión 3/4 motores).

Guía XFI completa: **[LEEME_XFI.md](./LEEME_XFI.md)**  
Acceso rápido monorepo: **[../../XFI.md](../../XFI.md)**

## Qué hay

| Archivo | Rol |
|---------|-----|
| `01_CONTEXT_ZYPYZAPE.txt` | Parámetros twin eólico + Quijote |
| `02_MATH_QUIJOTE_3vs7.txt` | 3 pales vs 7, J constante |
| `03_TASCA_GPT.txt` | Tareas de validación |
| `INICIO_avion_4turbinas_extract.txt` | Idea original chat (avión/coche) |
| `gemelo_xfi_avion_4turbinas.py` | **Gemelo XFI** (default N=3) |
| `zypy_zape_digital_twin_INICIO.py` | Primer twin N=5 |
| `zypyzape_twin_v2.py` | Twin documentado 5 nodos |
| `xfi_avion_dinamica_N*.png` | Gráficas de simulación |

## XFI — ejecutar

```bat
python gemelo_xfi_avion_4turbinas.py
python gemelo_xfi_avion_4turbinas.py --n 3
python gemelo_xfi_avion_4turbinas.py --compare
```

| N | Arquitectura |
|---|--------------|
| **3** | 1 gen + 1 thr + 1 buf (120°) — avión del chat |
| **4** | rueda gen/thr/buf — más potencia |

- **Quijote** = masas en palas (“molino en molino”).  
- **XFI** = control captura↔empuje; **no** perpetuum (η&lt;1).  
- Parte del **1** del trato **33×1**.
