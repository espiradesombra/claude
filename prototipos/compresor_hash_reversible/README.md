# Compresor por enlaces de hash reversibles (VHC2)

Hay **varios** hashes reversibles en el monorepo; este prototipo los **anida**.

Ver [`MAPA_HASHES_REVERSIBLES.md`](MAPA_HASHES_REVERSIBLES.md).

## Pipeline por defecto

```text
bloque → cinematic → espejo → toffoli → k3 → outer
         (meta)      (R,θ,t)   (puerta)  (huella)
```

```bat
cd prototipos\compresor_hash_reversible
python compresor_enlaces.py families
python compresor_enlaces.py pack ..\..\README.md -o demo.vhc -b 1024
python compresor_enlaces.py pack archivo.bin -o out.vhc -p cinematic,espejo
python compresor_enlaces.py verify demo.vhc
python compresor_enlaces.py unpack demo.vhc -o README.out.md
```

Contenedor `.vhc` = ZIP (`manifest.json` + `dict/<sha256>`). Formato **VHC2**.

## Siguiente

- Capa `oceanico` + `wifi` (params físicos).
- Modo generativo (sin diccionario) cuando el bloque sea regenerable solo con meta.
- Plugin AntiPC `HASH_CHAIN_PACK`.
