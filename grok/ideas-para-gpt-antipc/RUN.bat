@echo off
chcp 65001 >nul
cd /d "%~dp0"
echo.
echo  ╔══════════════════════════════════════╗
echo  ║  VMA JUST RUN — toolkit unificado    ║
echo  ╚══════════════════════════════════════╝
echo.
if exist "bin\vma-run.exe" (
  if "%~1"=="" (
    bin\vma-run.exe demo
  ) else (
    bin\vma-run.exe %*
  )
) else (
  if "%~1"=="" (
    python vma-run\vma_run.py demo
  ) else (
    python vma-run\vma_run.py %*
  )
)
echo.
pause