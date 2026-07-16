@echo off
cd /d "%~dp0.."
echo === AntiPC CMD v0.2 ===
echo.
antipc.cmd version
echo.
antipc.cmd hash --text "Hola VMA"
echo.
antipc.cmd mdc factor 10403
echo.
antipc.cmd wave --samples 3
echo.
antipc.cmd mechanical 2 3
echo.
pause