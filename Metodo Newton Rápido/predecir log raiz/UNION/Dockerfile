FROM ubuntu:22.04

# Herramientas de build y libs
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential cmake git python3 python3-pip \
    libmpfr-dev libgmp-dev && \
    rm -rf /var/lib/apt/lists/*

# Python para generar informe
RUN pip3 install numpy scipy matplotlib pandas jinja2

# Descargar fuentes (SLEEF, CRlibm)
WORKDIR /src
RUN git clone --depth=1 https://github.com/shibatch/sleef.git
RUN git clone --depth=1 https://github.com/SixTrack/crlibm.git

# Construir CRlibm
WORKDIR /src/crlibm
RUN ./configure --prefix=/usr && make -j && make install

# Construir SLEEF
WORKDIR /src/sleef/build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DSLEEF_SHOW_CONFIG=ON -DBUILD_SHARED_LIBS=ON .. && \
    make -j && make install

# Copiar harness y scripts
WORKDIR /work
COPY CMakeLists.txt /work/CMakeLists.txt
COPY bench_log.cpp   /work/bench_log.cpp
COPY bench_root.cpp  /work/bench_root.cpp
COPY tu_algo.hpp     /work/tu_algo.hpp
COPY tu_algo.cpp     /work/tu_algo.cpp
COPY gen_report.py   /work/gen_report.py
COPY README.md       /work/README.md

# Compilar y ejecutar
RUN cmake -B build -S . -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
CMD ["/bin/bash","-lc","./build/bench_log && ./build/bench_root && python3 gen_report.py && mkdir -p /out && cp -a report.html results/*.csv /out/"]
