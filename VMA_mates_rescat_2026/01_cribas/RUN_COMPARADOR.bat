@echo off
chcp 65001 >nul
title Comparador de cribas VMA
cd /d "%~dp0"

echo.
echo  Comparador de cribas VMA + Eratostenes
echo  Teoria: TEORIA_CRIBAS.md
echo.
python benchmark_cribas.py --quick
if errorlevel 1 (
  echo.
  echo  [AVISO] Alguna criba no coincidio con la referencia.
)
echo.
echo  Resultados en carpeta results\
echo.
pause
