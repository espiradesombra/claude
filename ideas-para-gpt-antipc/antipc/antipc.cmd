@echo off
chcp 65001 >nul
if exist "%~dp0dist\antipc.exe" (
    echo %*| findstr /i /c:"mdc analyze" /c:"mdc visual" /c:"mk hash" /c:"network demo" >nul
    if not errorlevel 1 goto :pycli
    if exist "%~dp0dist\k3hash.dll" set "PATH=%~dp0dist;%PATH%"
    "%~dp0dist\antipc.exe" %*
    exit /b %ERRORLEVEL%
)
:pycli
cd /d "%~dp0src\antipc"
python cli.py %*
exit /b %ERRORLEVEL%