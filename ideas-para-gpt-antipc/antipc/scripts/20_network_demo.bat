@echo off
chcp 65001 >nul
cd /d "%~dp0..\src\antipc"
echo AntiPC network demo — B slot-ring ^> D Grafcet
python cli.py network demo --hubs 4 --duration 2.5 --packets 15000
exit /b %ERRORLEVEL%