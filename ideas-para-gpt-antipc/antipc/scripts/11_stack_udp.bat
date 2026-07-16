@echo off
chcp 65001 >nul
title AntiPC — Stack + UDP HMAC
cd /d "%~dp0..\src\antipc"

echo.
echo  AntiPC — AntiPCStack enlazado a hubs UDP (HMAC HASHTOOLCODE)
echo.
python demo_stack_udp.py --root "%~dp0.." --hubs 4 --max-files 40
echo.
pause