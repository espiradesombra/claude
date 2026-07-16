# Aleatorovix + Firma Geométrica

Hub de documentación para el motor de **convergencia geométrica** (pares `10/11/00/01` → perímetro) y el organismo **Aleatorovix** (entropía + máscara Lila + cifrado por fase π/e).

**Autor:** Víctor Manzanares Alberola · **Proyecto:** VMA / 33×1 · **Uso civil únicamente**

Índice detallado: [`INDICE_ALEATOROVIX_GEO.txt`](INDICE_ALEATOROVIX_GEO.txt)

---

## Rutas canónicas

| Componente | Ruta en repo |
|------------|--------------|
| Teoría + demos Python | `encriptacionGeometrica/` |
| Organismo Aleatorovix | `aleatorovix/` |
| Núcleos C (copia auditable) | `docs/aleatorovix-firma-geometrica/fuentes-c/` |
| DLL / CLI integrada | `antipc/src/antipc/` |
| Teorema fase K3 | `teoremas/19_teorema_phaseamplifier_k3.txt` |

---

## Comandos AntiPC (CLI v0.14)

Desde `antipc/src/antipc`:

```bat
python cli.py version

REM Convergencia geométrica binaria (C, requiere DLL)
python cli.py geo --bits 0101101011000010 --tales 3,5,8,13,21 --puntos 6,12,18

REM Aleatorovix + fase geométrica masiva (round-trip demo)
python cli.py geo-masivo --text "hola vma" --semilla 43210

REM Motor acordeón K3 XOR (proyecto 33×1)
python cli.py k3 stream-xor --text "libro4-33x1" --base 33 --rel 1

REM Vía corpus Libro 4
python cli.py libro run L4-convergencia
python cli.py libro run L4-aleatorovix
python cli.py libro run L4-k3-stream

python cli.py teorema info 19
```

Compilar DLL nativa (Windows):

```bat
cd antipc\scripts
21_build_antipc_native.bat
```

---

## Demos Python (sin DLL)

```bat
cd encriptacionGeometrica
python _demo_convergencia.py

cd ..\aleatorovix
python aleatorovix.py -n 10
python demos\demo_organismo.py
python benchmarks\benchmark_aleatorovix.py
```

---

## Modelo operativo (resumen)

### Firma geométrica

1. Estado compartido: secuencia Thales, tipos de figura, puntos base.
2. Lectura de pares de bits; huérfano final → se añade `0`.
3. Cada par aporta lados al perímetro según la tabla `10/11/00/01`.
4. Salida: `Decimal` de alta precisión (perímetro de colapso). Reversible con el mismo estado.

### Aleatorovix (Gemini / port C)

1. Semilla 16 bits → cifrado vía perímetro × π ofuscado × e convergente.
2. Chorro masivo de bits; XOR con datos.
3. Descifrado: backtracking sobre residuo decimal (Python) o hash de fase (C).

### K3 industrial (`k3_core.c`)

- Motor acordeón: `L += v+1`, `v += 2`; validación `5*L <= 2*v+1`.
- Stream: `(L ^ v) * 0x9E3779B97F4A7C15` → XOR byte a byte.
- Parámetros 33×1: `base=33`, `rel=1`.

---

## Advertencias

- `k3_launcher.py` **elimina el archivo original** tras cifrar. Usar solo copias de prueba.
- `k3_core.c` usa `mlock`/`memset` orientado a Linux; en Windows compilar en WSL o entorno de prueba.
- Motor toy / auditable — no sustituto de criptografía estándar (AES, etc.).

---

## Enlaces

- Repo: https://github.com/espiradesombra/claude
- Integración C: [`docs/integracion-antipc/`](../integracion-antipc/)
- Mapa port: [`01_MAPA_PORT_A_C.md`](../integracion-antipc/01_MAPA_PORT_A_C.md)