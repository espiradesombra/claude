@echo off
chcp 65001 >nul
title AntiPC — Benchmark A-D
cd /d "%~dp0..\src\antipc"

echo.
echo  AntiPC — Benchmark arquitecturas A, B, C, D
echo.
python benchmark.py --packets 2000 --payload 128 --hubs 4
echo.
pause