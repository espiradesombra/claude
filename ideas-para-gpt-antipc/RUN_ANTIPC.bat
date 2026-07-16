@echo off
chcp 65001 >nul
cd /d "%~dp0"
echo === AntiPC — benchmark UDP 4 hubs (localhost) ===
if exist "bin\vma-run.exe" (
  bin\vma-run.exe antipc
) else (
  python vma-run\vma_run.py antipc
)
pause