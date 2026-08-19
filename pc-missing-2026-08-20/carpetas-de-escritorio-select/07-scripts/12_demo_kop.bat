@echo off
chcp 65001 >nul
cd /d "%~dp0..\src\antipc"
echo === AntiPC KOP binario (K3 MicroKernel) ===
python demo_kop_binary.py
pause