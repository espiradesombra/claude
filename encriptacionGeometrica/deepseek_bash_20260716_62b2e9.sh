#!/bin/bash
mkdir -p build
gcc -O3 -march=native -fPIC -shared -o build/libk3.so src/k3.c -lm
cp python/k3.py build/
cp python/demo.py build/
echo "✅ Build complete. Ejecuta: cd build && python3 demo.py"