@echo off
chcp 65001 >nul
cd /d "%~dp0.."
echo === MDC visual — dos trenes ===
echo.
set N=10403
if not "%~1"=="" set N=%~1
echo n=%N%
echo.
antipc.cmd mdc analyze %N% --proper
echo.
antipc.cmd mdc visual %N% --proper -o "output\mdc_trains_%N%.gif"
echo.
echo GIF en output\mdc_trains_%N%.gif
echo Ventana animada: antipc.cmd mdc visual %N% --proper --gui
pause