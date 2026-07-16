@echo off
chcp 65001 >nul
title AntiPC — Arquitectura Unificada
cd /d "%~dp0..\src\antipc"

echo.
echo  AntiPC — Arquitectura HASHTOOLCODE + MDC + Flow Kernel + Red
echo.
python demo_arquitectura.py --root "%~dp0.."
echo.
pause