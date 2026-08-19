@echo off
chcp 65001 >nul
title AntiPC v0.1.0 — Flow Kernel
cd /d "%~dp0..\src\antipc"

echo.
echo  AntiPC v0.1.0 — Flow Kernel (gptcomputing.txt roadmap)
echo.
python demo_v010.py --mode battery --count 500 --repeat 0.35
echo.
pause