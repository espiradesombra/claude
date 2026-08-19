@echo off
cd /d "%~dp0.."
python benchmarks\benchmark_aleatorovix.py
if errorlevel 1 exit /b 1
echo.
type benchmarks\RESULTADOS_BENCHMARK.txt
pause