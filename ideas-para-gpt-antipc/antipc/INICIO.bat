@echo off
chcp 65001 >nul
title AntiPC Sale-It v0.1
cd /d "%~dp0"

:menu
cls
echo.
echo   ╔══════════════════════════════════════════════════════════╗
echo   ║  AntiPC — Adaptive Network Through Parallel Computing   ║
echo   ║  Sale-It v0.1                                           ║
echo   ╚══════════════════════════════════════════════════════════╝
echo.
echo   [1] Instalar / verificar Python
echo   [2] Benchmark A-D (local)
echo   [3] Benchmark A-E con UDP real  ^(recomendado para demo^)
echo   [4] Solo UDP ^(4 hubs^)
echo   [5] Verificar motor K3
echo   [6] Abrir carpeta del producto
echo   [7] Leer LEEME.txt
echo   [8] Leer PRODUCT_SHEET_EN.txt
echo   [9] Demo v0.1.0 Flow Kernel ^(gptcomputing.txt^)
echo   [I] Demo INDUSTRIAL ^(auditoria ficheros + K(N)^)
echo   [A] Arquitectura UNIFICADA ^(HASHTOOLCODE+MDC+AntiPC^)
echo   [U] Stack + UDP enlazado ^(HMAC + hubs remotos^)
echo   [0] Salir
echo.
set /p opcion="  Elige opcion: "

if "%opcion%"=="1" call scripts\01_instalar.bat & goto menu
if "%opcion%"=="2" call scripts\02_benchmark.bat & goto menu
if "%opcion%"=="3" call scripts\03_benchmark_udp.bat & goto menu
if "%opcion%"=="4" call scripts\04_udp_solo.bat & goto menu
if "%opcion%"=="5" call scripts\05_verificar_k3.bat & goto menu
if "%opcion%"=="6" call scripts\06_abrir_carpeta.bat & goto menu
if "%opcion%"=="7" start notepad LEEME.txt & goto menu
if "%opcion%"=="8" start notepad PRODUCT_SHEET_EN.txt & goto menu
if "%opcion%"=="9" call scripts\07_demo_v010.bat & goto menu
if /i "%opcion%"=="I" call scripts\09_demo_industrial.bat & goto menu
if /i "%opcion%"=="A" call scripts\10_arquitectura.bat & goto menu
if /i "%opcion%"=="U" call scripts\11_stack_udp.bat & goto menu
if "%opcion%"=="0" exit /b 0
goto menu