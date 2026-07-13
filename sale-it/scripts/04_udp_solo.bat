@echo off
chcp 65001 >nul
title AntiPC — UDP solo
cd /d "%~dp0..\src\antipc"

echo.
echo  AntiPC — Arquitectura E (UDP real, 4 hubs)
echo.
python udp_benchmark.py --packets 2000 --payload 64 --hubs 4
echo.
pause