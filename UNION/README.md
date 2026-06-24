# logbench: Comparativa de tu algoritmo vs TOP libm

Este paquete ejecuta benchmarks de **logaritmos** y **raíces** comparando:

- **TU_ALGO** (doble ajuste simultáneo)
- **libm** del sistema (`std::log`, `std::pow`)
- **CRlibm** (log natural con *correct rounding*)
- **SLEEF** (≈1 ULP, vectorizable)

Y genera `report.html` con **tiempo** y **ULP@95%**.

## Docker (recomendado)
```bash
docker build -t logbench .
docker run --rm -it -v "$PWD/out:/out" logbench
# Ver ./out/report.html
```

## Instalación nativa (Ubuntu/Debian)
```bash
sudo apt-get update && sudo apt-get install -y build-essential cmake git libmpfr-dev libgmp-dev python3 python3-pip
pip3 install numpy scipy matplotlib pandas jinja2

git clone --depth=1 https://github.com/SixTrack/crlibm.git && (cd crlibm && ./configure && make -j && sudo make install)
git clone --depth=1 https://github.com/shibatch/sleef.git && (cd sleef && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j && sudo make install)

cmake -B build -S . -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/bench_log && ./build/bench_root && python3 gen_report.py
# Ver report.html
```

## Notas
- CRlibm: log con *correct rounding* en doble precisión.
- SLEEF: variantes 1-ULP con rutas SIMD.
- El oráculo de alta precisión con MPFR no está incluido en C++ (para portabilidad); `ULP` se calcula contra `libm`/CRlibm. Para una auditoría formal, amplíe con MPFR.
