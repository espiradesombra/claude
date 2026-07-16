@echo off
setlocal
cd /d "%~dp0\.."
echo.
echo  AntiPC — subida a GitHub desde Escritorio
echo  Carpeta: %CD%
echo.

where git >nul 2>&1 || (echo ERROR: instala Git for Windows & exit /b 1)

if not exist ".git" (
  echo [1/5] Inicializando git...
  git init
  git branch -M main
) else (
  echo [1/5] Repo git ya existe.
)

echo [2/5] Anadiendo archivos...
git add -A
git status -sb

echo.
set /p MSG= Mensaje commit [AntiPC v0.14]: 
if "%MSG%"=="" set MSG=AntiPC v0.14 — DeepSeek, jerk, desmemoriada, teoremas CLI

git commit -m "%MSG%" 2>nul
if errorlevel 1 (
  echo Sin cambios nuevos o commit ya hecho.
)

git remote get-url origin >nul 2>&1
if errorlevel 1 (
  echo.
  echo Remote no configurado. Ejemplos:
  echo   https://github.com/espiradesombra/antipc.git
  echo   https://github.com/espiradesombra/claude.git
  set /p REMOTE= URL origin: 
  if not "!REMOTE!"=="" git remote add origin "!REMOTE!"
)

echo [3/5] Push a origin main...
echo      Usa tu PAT de GitHub como contraseña si lo pide.
git push -u origin main

if errorlevel 1 (
  echo.
  echo PUSH FALLO — revisa SUBIDA_GITHUB.txt
  echo  - Token PAT con permiso repo
  echo  - Repo creado en github.com/new
  exit /b 1
)

echo.
echo OK — subido a GitHub.
endlocal