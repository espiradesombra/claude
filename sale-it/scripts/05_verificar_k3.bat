@echo off
chcp 65001 >nul
title AntiPC — Verificar K3
cd /d "%~dp0..\src\antipc"

echo.
echo  AntiPC — Verificacion motor K3
echo.
python verify_k3.py
echo.
pause