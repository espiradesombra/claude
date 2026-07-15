@echo off
chcp 65001 >nul
cd /d "%~dp0"
set OUT=%~dp0just-run-unificado.zip
set STAGE=%TEMP%\just-run-stage-%RANDOM%

echo === EMPACAR just run unificado ===
echo.

if exist "%OUT%" del /f /q "%OUT%"
if exist "%STAGE%" rmdir /s /q "%STAGE%"
mkdir "%STAGE%"

echo [1/3] Copiando paquete (excluyendo cache y build)...
robocopy "%~dp0" "%STAGE%\just-run" /E /XD __pycache__ build vma_k3.egg-info /XF *.pyc *.pyo "Nueva carpeta comprimida (en zip).zip" /NFL /NDL /NJH /NJS /nc /ns /np >nul

echo [2/3] Comprimiendo...
powershell -NoProfile -Command "Compress-Archive -Path '%STAGE%\just-run\*' -DestinationPath '%OUT%' -Force"

echo [3/3] Limpiando temporal...
rmdir /s /q "%STAGE%"

for %%A in ("%OUT%") do set SIZE=%%~zA
set /a SIZE_MB=%SIZE%/1048576
echo.
echo OK: %OUT%
echo Tamano aprox: %SIZE_MB% MB (sin exe grande puede ser 3-8 MB)
echo.
echo Incluye:
echo   - README.md, VALOR_ZIP_UNIFICADO.txt, TABLA_CONTENIDO_VMA.md
echo   - gptcomputing\CONTINUACION_PARA_GPT.txt  ^(para ChatGPT^)
echo   - antipc\, vma-run\, gptcomputing\, archivos-vma\, ...
echo.
echo Para ChatGPT: adjunta zip o solo gptcomputing\ + antipc\src\ + docs\
pause