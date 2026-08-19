# k3hash — Librería de hash no criptográfico (Windows / macOS / Linux)

Basada en tu algoritmo K3 original (mezcla por rotaciones + XOR + constante
de dispersión `0x9E3779B1`). Empaquetada como librería C portable con CMake,
para compilar nativamente en los tres sistemas sin tocar el código.

**Importante:** esto es un hash de propósito general (32 bits), no
criptográfico. Sirve para checksums de integridad, deduplicación, tablas
hash, fingerprinting y sharding — no para contraseñas ni autenticación
(para eso usa la librería HMAC-SHA256 ya entregada).

## Estructura
```
k3hash/
├── include/k3hash.h      ← API pública
├── src/k3hash.c          ← implementación (algoritmo original + streaming)
├── examples/k3hash_cli.c ← CLI de ejemplo
└── CMakeLists.txt        ← build multiplataforma
```

## Herramientas incluidas (los tres portables)

Todas se compilan automáticamente al hacer `cmake --build build` — no hace
falta nada extra.

### 1. `k3dedup` — duplicados exactos
```
find /mi/carpeta -type f | ./k3dedup
```
Agrupa por (tamaño + hash K3). Solo detecta copia byte a byte idéntica.

### 2. `k3simil` — similitud / plagio parcial
```
find /mi/carpeta -name "*.txt" | ./k3simil 0.30
```
Trocea cada documento en frases solapadas de 5 palabras ("shingles"),
hashea cada una con K3, y compara los conjuntos con el índice de Jaccard.
El argumento opcional es el umbral (0.0–1.0, por defecto 0.30). Detecta
reformulaciones parciales, no solo copia literal.

### 3. `k3search` — buscador de contenido
```
find /mi/biblioteca -name "*.txt" | ./k3search index biblioteca.k3idx
./k3search query biblioteca.k3idx palabra1 palabra2
```
Construye un índice invertido (término → documentos), usando K3 para
ubicar cada término en una tabla hash O(1). El índice se guarda en texto
plano (`.k3idx`), portable entre los tres sistemas operativos sin
conversión. La consulta puntúa documentos por frecuencia de los términos
buscados y muestra los 20 mejores.

En Windows, sustituye `find ... -type f` por `dir /s /b <carpeta>`.

## API principal

```c
K3HashConfig cfg = k3_config_default();

/* Un solo bloque de memoria */
uint32_t h = k3_hash_buffer(datos, longitud, &cfg);

/* Streaming (archivos grandes, datos por partes) */
K3HashCtx ctx;
k3_hash_init(&ctx, &cfg);
k3_hash_update(&ctx, parte1, len1);
k3_hash_update(&ctx, parte2, len2);
uint32_t h = k3_hash_final(&ctx);

/* Fichero directo */
uint32_t h;
k3_hash_file("ruta/al/archivo", &cfg, &h);
```

---

## Compilar en Windows

### Opción A: Visual Studio (recomendado)
1. Instala **CMake** (cmake.org) y **Visual Studio 2022** con "Desktop
   development with C++".
2. Desde "Developer Command Prompt for VS":
   ```
   cmake -S . -B build -G "Visual Studio 17 2022" -A x64
   cmake --build build --config Release
   ```
   Esto genera `build\Release\k3hash.dll`, `k3hash.lib` (import lib) y
   `k3hash_cli.exe`.
3. Para usar la DLL en tu propio proyecto: incluye `k3hash.h`, define
   `K3HASH_SHARED` antes de incluirlo, enlaza `k3hash.lib` y copia
   `k3hash.dll` junto a tu `.exe`.

### Opción B: MinGW (sin Visual Studio)
```
cmake -S . -B build -G "MinGW Makefiles"
cmake --build build
```

---

## Compilar en macOS
```
cmake -S . -B build
cmake --build build
```
Esto genera `build/libk3hash.dylib` (o `libk3hash.a` con
`-DK3HASH_BUILD_SHARED=OFF`) y `build/k3hash_cli`.

Para enlazarlo en Xcode: añade `k3hash.dylib` a "Link Binary With
Libraries" y la carpeta `include/` a "Header Search Paths".

---

## Compilar en Linux (ya probado en este entorno)
```
cmake -S . -B build
cmake --build build
./build/k3hash_cli --text "hola"
./build/k3hash_cli --file /ruta/a/archivo
```

## Instalar en el sistema (los tres SO)
```
cmake --install build --prefix /ruta/de/instalacion
```

## Nota sobre reproducibilidad
El hash es determinista: mismo input + misma `K3HashConfig` → mismo
resultado, siempre, en Windows/macOS/Linux (no depende de endianness de la
máquina porque los bloques se ensamblan explícitamente byte a byte).
