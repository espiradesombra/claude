@echo off
cd /d "%~dp0"
echo === ALEATOROVIX v1.0 ===
echo.
echo [1] Organismo 10 ciclos
echo [2] Demo primos + criba
echo [3] Demo MDC N=10403
echo [4] Organismo 30 ciclos con MDC
echo [5] Benchmark completo
echo.
set /p OPCION=Elige 1-5 (Enter=1): 
if "%OPCION%"=="" set OPCION=1
if "%OPCION%"=="1" python aleatorovix.py
if "%OPCION%"=="2" python demos\demo_primos.py
if "%OPCION%"=="3" python demos\demo_mdc_N.py
if "%OPCION%"=="4" python aleatorovix.py -n 30 --n-mdc 10403
if "%OPCION%"=="5" python benchmarks\benchmark_aleatorovix.py
echo.
pause