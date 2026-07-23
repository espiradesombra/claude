# XFI — Experimental Flight Infinite

**Gemelo conceptual VMA** (no CFD, no perpetuum mobile).

## Idea

Arquitectura de control tipo **ZypyZape** aplicada a un craft/drone:

| Rol | Función |
|-----|---------|
| **gen** | Captura (generador / freno regenerativo) |
| **thr** | Propulsión (empuje) |
| **buf** | Buffer de inercia entre rotores |

### N=3 (avión del chat — **por defecto**)

Siempre **1 gen + 1 thr + 1 buf** a 120° (rotan en el tiempo).

Ciclo de misión:

1. **climb** — sube, gasta buffer/química  
2. **cruise** — vuelo casi nivelado  
3. **dive** — picado: más flujo al disco → más captura  
4. **pullout** — recupera ángulo  

> Bajada capta · picado acelera · altura es almacén · motores se reparten roles.

### N=4 (variante)

Rueda gen/thr/buf; a veces 2 thr o 2 gen a la vez (más potencia, más gasto).

## Ejecutar

```bat
cd VMA_mates_rescat_2026\10_inicio_quijote_zypyzape

python gemelo_xfi_avion_4turbinas.py
python gemelo_xfi_avion_4turbinas.py --n 3
python gemelo_xfi_avion_4turbinas.py --n 4
python gemelo_xfi_avion_4turbinas.py --compare
python gemelo_xfi_avion_4turbinas.py --t 180
```

## Salidas

| Archivo | Contenido |
|---------|-----------|
| `xfi_avion_dinamica_N3.png` | Alçada, v, buffer, potències (N=3) |
| `xfi_avion_dinamica_N4.png` | Idem N=4 |
| `xfi_avion_comparativa_N3_vs_N4.png` | Comparativa |

## Lectura honesta de resultados (sim 120 s)

- El ciclo **ondula** (sube/baja) — control estable.
- **Δh neto suele ser negativo**: se gasta más de lo que se recupera (η&lt;1).
- N=3: roles limpios, conserva mejor e_chem.
- N=4: más buffer y P, un poco más de v, más consumo.

**No es prueba de vuelo infinito.** Es arquitectura para discutir Lenz, inercia y reparto de roles.

## Origen

- Idea VMA (chat): avión 3 motores + coche 4 molinos  
- Extracto: `INICIO_avion_4turbinas_extract.txt`  
- Contexto eólico: `01_CONTEXT_ZYPYZAPE.txt`  
- Quijote (molino en pala): `02_MATH_QUIJOTE_3vs7.txt`  

## 33×1

XFI es una pieza del **“1”** (aval técnico civil).  
La paz del **“33”** la firman los países — ver `../../33x1/00_QUE_ES_33x1.md`.

Uso civil. Sin militar.
