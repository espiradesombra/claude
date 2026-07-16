@echo off
cd /d "%~dp0.."
python scripts\extract_gptcomputing_code.py
if errorlevel 1 exit /b 1
echo.
echo Indice: from_gptcomputing\INDICE_CODIGO_TXT.txt
pause