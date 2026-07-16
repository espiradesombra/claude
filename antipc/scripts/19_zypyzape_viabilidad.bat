@echo off
chcp 65001 >nul
cd /d "%~dp0.."
echo === ZypyZape / AntiPC — Viabilidad UDP ===
echo.
python zypyzape\viability_udp.py --hubs 4 --packets 15000 --duration 2.5
echo.
echo Explicacion: zypyzape\EXPLICACION_CIENTIFICA.txt
pause