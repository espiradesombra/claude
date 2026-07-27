# Matemáticas (eje)

Teoría de números y métodos VMA: cribas, MDC, Newton Rápido, MRAUV, Goldbach, Sophie Germain, fichas de teoremas.

## Carpetas clave

| Carpeta | Contenido |
|---------|-----------|
| `VMA_mates_rescat_2026/` | Pack limpio 2026-07-23: 00 sunraman … 10 XFI (cada subcarpeta con README) |
| `vma-methods/` | Librería ejecutable (cmd + paquete Python) |
| `teoremas/` | Fichas 01–30+, formal / ingeniería / espacial |
| `teoremasgrok/` | Línea Grok de las fichas |
| `Metodo Newton Rápido/` | Papers, oráculos MEcuation, benches |
| `Metodo densidad MRAUV/`, `mrauv/` | Densidad de primos |
| `Goldbach/`, `Sophie Germain/`, `siguiente_primo/` | Líneas clásicas |
| `COMPARATIVA/`, `graficas y explicaciones/` | Comparativas y figuras comentadas |
| `Libro5 Factorizacion con 2v+3/` | Corpus L5 / MDC |

## Pruebas rápidas (pack rescat)

```bat
cd VMA_mates_rescat_2026\01_cribas
python -c "from cribas import comparar_cribas; comparar_cribas(500)"

cd ..\05_sofi_fermat_goldbach
python -c "from siguiente_primo import validate; validate(40)"
```

## Lectura

- README de cada carpeta  
- Fichas `teoremas/01_…`  
- Wiki [Glosario](Glosario) (MDC, criba, MRAUV…)
