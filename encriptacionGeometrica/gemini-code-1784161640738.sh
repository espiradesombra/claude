# Esquema de compilación y ejecución industrial consolidado:
gcc -O3 -march=native -shared -fPIC -o nucleo_k3.so nucleo_k3.c -lm
python3 aleatorovix_final.py