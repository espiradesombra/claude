@echo off
chcp 65001 >nul
cd /d "%~dp0aleatorovix"
echo === ALEATOROVIX — just run ===
echo.
echo [1] Criba Desmemoriada (Python)
python criba_desmemoriada.py
echo.
echo [2] MDC v6 memoria adaptativa (demo)
python mdc_v6_aleatorovix.py
echo.
echo [3] Organismo Lila C — requiere WSL/Linux:
echo     wsl gcc -O2 -o aleatorovix organismo_lila_v99.c -lm ^&^& ./aleatorovix
echo.
echo Lee CONCEPTO.txt y ALEATOROVIX.txt para contexto completo.
pause