@echo off
chcp 65001 >nul
cd /d "%~dp0..\src\antipc"
echo AntiPC — inventario industrial
python -m cli industrial inventory
python -m cli industrial audit --max-files 50 --full
exit /b %ERRORLEVEL%