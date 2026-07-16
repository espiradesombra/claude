@echo off
chcp 65001 >nul
title Compilar k3hash.dll — HASHTOOLCODE
cd /d "%~dp0..\native\k3hash"

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
    echo Instala Visual Studio con "Desktop development with C++" o CMake standalone.
    exit /b 1
)

echo CMake: %CMAKE%
echo.

if exist build rmdir /s /q build

"%CMAKE%" -S . -B build -G "Visual Studio 18 2026" -A x64 -DK3HASH_BUILD_EXAMPLES=OFF
if errorlevel 1 (
    echo Intentando Visual Studio 17 2022...
    if exist build rmdir /s /q build
    "%CMAKE%" -S . -B build -G "Visual Studio 17 2022" -A x64 -DK3HASH_BUILD_EXAMPLES=OFF
)
if errorlevel 1 (
    echo Intentando generador por defecto...
    if exist build rmdir /s /q build
    "%CMAKE%" -S . -B build -DK3HASH_BUILD_EXAMPLES=OFF
)

"%CMAKE%" --build build --config Release
if errorlevel 1 (
    echo [ERROR] Compilacion fallida
    exit /b 1
)

set "DLL=build\Release\k3hash.dll"
if not exist "%DLL%" set "DLL=build\k3hash.dll"
if not exist "%DLL%" (
    echo [ERROR] k3hash.dll no generada
    exit /b 1
)

copy /Y "%DLL%" "..\..\k3hash.dll" >nul
copy /Y "%DLL%" "..\..\src\antipc\k3hash.dll" >nul
copy /Y "%DLL%" "..\..\src\antipc\native\k3hash\build\Release\k3hash.dll" >nul
if exist "..\..\dist" copy /Y "%DLL%" "..\..\dist\k3hash.dll" >nul

echo.
echo [OK] k3hash.dll desplegada en:
echo   antipc\k3hash.dll
echo   src\antipc\k3hash.dll
echo   src\antipc\native\k3hash\build\Release\
echo   dist\  (si existe)
echo.
cd /d "%~dp0..\src\antipc"
python verify_k3.py
exit /b %ERRORLEVEL%