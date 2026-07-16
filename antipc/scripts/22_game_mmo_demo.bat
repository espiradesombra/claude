@echo off
chcp 65001 >nul
cd /d "%~dp0..\src\antipc"
echo AntiPC MMO — B slot-ring ^> D Grafcet ^> WorldState
python -m cli game demo --players 128 --duration 5 --shards 4
exit /b %ERRORLEVEL%