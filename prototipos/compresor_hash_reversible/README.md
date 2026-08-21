# Compresor por enlaces de hash reversibles (prototipo)

Parte de:

- `hash cinematico de otras formas y calculadoras/` (1 acelera, 0 decel+racha, metadatos)
- `compresor y ultima calculadora/` (hash + params para reconstruir)
- K3 (`hashtool-work` / AntiPC)

## Uso

```bat
cd prototipos\compresor_hash_reversible
python compresor_enlaces.py pack ..\..\README.md -o demo.vhc -b 1024
python compresor_enlaces.py verify demo.vhc
python compresor_enlaces.py unpack demo.vhc -o README.out.md
```

Contenedor `.vhc` = ZIP (`manifest.json` + `dict/<sha256>`).

## Qué hace

1. Trocea el fichero.
2. Por bloque: hash **cinemático** + **metadatos** + digest **K3** encadenado  
   `H_i = K3(H_{i-1} || bloque || digest_cin)`.
3. Deduplica bloques idénticos (mismo SHA-256) → solo guarda un payload.
4. Unpack/verify **recalcula** cinemática + cadena (reversible en el sentido de verificación+reproducción de la trayectoria).

## Qué no es (aún)

- No regenera bytes solo desde el digest sin diccionario (haría falta un modo “espejo generativo”).
- No sustituye zstd en datos aleatorios.
- La calculadora WiFi / AntiPC se integra en el **siguiente** prototipo.

## Siguiente

- Modo espejo: bloques generables solo con `meta` (sin store).
- Enchufar al runtime AntiPC como plugin `HASH_CHAIN_PACK`.
- Calculadora WiFi/ondas como fuente de entropía / params de pasada.
