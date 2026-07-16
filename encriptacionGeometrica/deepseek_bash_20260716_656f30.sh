#!/bin/bash
# Compila la biblioteca K3 en C y la instala en ./build/

mkdir -p build
gcc -O3 -march=native -fPIC -shared -o build/libk3.so src/k3_engine.c -lm
echo "✅ Biblioteca compilada en build/libk3.so"

# Copiar cabecera
cp include/k3_engine.h build/

# Crear enlace simbólico para Python
ln -sf $(pwd)/build/libk3.so python/libk3.so
echo "✅ Enlace creado para Python"