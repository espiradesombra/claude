@echo off
chcp 65001 >nul
title Compilar k3hash.dll — HASHTOOLCODE
cd /d "%~dp0"

where cmake >nul 2>&1
if errorlevel 1 (
    echo [ERROR] CMake no encontrado.
    echo Instala CMake + Visual Studio 2022 con "Desktop development with C++"
    echo Luego ejecuta desde "Developer Command Prompt for VS":
    echo   cmake -S . -B build -G "Visual Studio 17 2022" -A x64
    echo   cmake --build build --config Release
    pause
    exit /b 1
)

cmake -S . -B build -G "Visual Studio 17 2022" -A x64
if errorlevel 1 (
    echo Intentando generador por defecto...
    cmake -S . -B build
)
cmake --build build --config Release
if errorlevel 1 (
    echo [ERROR] Compilacion fallida
    pause
    exit /b 1
)

if exist "build\Release\k3hash.dll" (
    copy /Y "build\Release\k3hash.dll" "..\..\k3hash.dll"
    echo [OK] k3hash.dll copiada a antipc\
) else if exist "build\k3hash.dll" (
    copy /Y "build\k3hash.dll" "..\..\k3hash.dll"
    echo [OK] k3hash.dll copiada a antipc\
)

echo.
echo Ejecuta: python verify_k3.py
pause