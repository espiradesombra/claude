@echo off
chcp 65001 >nul
title AntiPC Sale-It — Instalacion
cd /d "%~dp0.."

echo.
echo  AntiPC Sale-It v0.1 — Verificacion de entorno
echo  =============================================
echo.

where python >nul 2>&1
if errorlevel 1 (
    echo  [ERROR] Python no encontrado. Instala Python 3.10+ desde python.org
    pause
    exit /b 1
)

python --version
echo.
python -c "import sys; assert sys.version_info>=(3,10), 'Se requiere Python 3.10+'; print('  [OK] Version compatible')"
if errorlevel 1 (
    echo  [ERROR] Version de Python insuficiente
    pause
    exit /b 1
)

echo.
echo  [OK] Entorno listo. Ejecuta 02_benchmark.bat o 03_benchmark_udp.bat
echo.
pause