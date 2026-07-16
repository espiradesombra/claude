@echo off
chcp 65001 >nul
title AntiPC — Demo Industrial
cd /d "%~dp0..\src\antipc"

echo.
echo  AntiPC — Auditoria industrial (K3_FILE + K3_DEDUP + Knowledge Resolver)
echo  Carpeta por defecto: sale-it (puedes pasar --root otra ruta)
echo.
python industrial_demo.py --root "%~dp0.." --max-files 80 --workers 4
echo.
pause