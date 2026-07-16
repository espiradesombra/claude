@echo off
chcp 65001 >nul
cd /d "%~dp0"
python -m vma_methods.cli %*
exit /b %ERRORLEVEL%