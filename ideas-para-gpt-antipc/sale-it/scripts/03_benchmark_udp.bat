@echo off
chcp 65001 >nul
title AntiPC — Benchmark A-E (UDP)
cd /d "%~dp0..\src\antipc"

echo.
echo  AntiPC — Benchmark completo A-E con UDP real
echo.
python benchmark.py --udp --packets 1000 --payload 64 --hubs 4
echo.
pause