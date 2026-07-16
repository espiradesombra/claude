@echo off
chcp 65001 >nul
title antipc-port-c — compilar prototipos (solo Escritorio)
cd /d "%~dp0"

set SRC=src\criba_modular6k.c src\criba_hibrida.c src\mdc_trains.c src\newton_rapido.c src\mdc_ksweep_predict.c src\port_demo.c
set INC=-Iinclude
set OUT=port_demo.exe

where cl >nul 2>&1
if %ERRORLEVEL%==0 (
    echo [MSVC] cl ...
    cl /nologo /O2 /W3 %INC% %SRC% /Fe:%OUT%
    if errorlevel 1 exit /b 1
    echo [OK] %OUT%
    %OUT%
    exit /b 0
)

where gcc >nul 2>&1
if %ERRORLEVEL%==0 (
    echo [GCC] gcc ...
    gcc -O3 -Wall %INC% %SRC% -o %OUT% -lm
    if errorlevel 1 exit /b 1
    echo [OK] %OUT%
    %OUT%
    exit /b 0
)

echo [ERROR] No hay cl ni gcc en PATH.
echo Abre "x64 Native Tools Command Prompt for VS" y ejecuta este .bat
exit /b 1