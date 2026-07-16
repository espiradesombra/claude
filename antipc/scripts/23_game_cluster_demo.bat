@echo off
chcp 65001 >nul
cd /d "%~dp0..\src\antipc"
echo AntiPC MMO CLUSTER — L3 hubs WORK/RESULT por shard
python -m cli game cluster --players 128 --shards 4 --hubs 4 --duration 5
exit /b %ERRORLEVEL%