@echo off
chcp 65001 >nul
title Compilar antipc_native.dll — núcleo C unificado
cd /d "%~dp0..\src\antipc\native\antipc_core"

set "CMAKE="
for %%P in (
    "C:\Program Files\Microsoft Visual Studio\18\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
    "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
    "C:\Program Files\CMake\bin\cmake.exe"
) do if exist %%P set "CMAKE=%%~P"

if not defined CMAKE (
    where cmake >nul 2>&1 && set "CMAKE=cmake"
)

if not defined CMAKE (
    echo [ERROR] CMake no encontrado.
    exit /b 1
)

echo CMake: %CMAKE%
if exist build rmdir /s /q build

"%CMAKE%" -S . -B build -G "Visual Studio 17 2022" -A x64
if errorlevel 1 (
    if exist build rmdir /s /q build
    "%CMAKE%" -S . -B build -DK3HASH_BUILD_EXAMPLES=OFF
)
if errorlevel 1 (
    echo [ERROR] Configuracion CMake fallida
    exit /b 1
)

"%CMAKE%" --build build --config Release
if errorlevel 1 (
    echo [ERROR] Compilacion fallida
    exit /b 1
)

set "DLL=build\Release\Release\antipc_native.dll"
if not exist "%DLL%" set "DLL=build\Release\antipc_native.dll"
if not exist "%DLL%" set "DLL=build\antipc_native.dll"
if not exist "%DLL%" (
    echo [ERROR] antipc_native.dll no generada
    exit /b 1
)

copy /Y "%DLL%" "..\..\antipc_native.dll" >nul
copy /Y "%DLL%" "..\..\..\antipc_native.dll" >nul
if exist "..\..\..\..\dist" copy /Y "%DLL%" "..\..\..\..\dist\antipc_native.dll" >nul

echo.
echo [OK] antipc_native.dll desplegada
echo   src\antipc\antipc_native.dll
echo   antipc\antipc_native.dll
echo.
cd /d "%~dp0..\src\antipc"
python -c "from native_engine import status_report; print(status_report())"
exit /b %ERRORLEVEL%