@echo off
chcp 65001 >nul
cd /d "%~dp0"
python -m vma_methods.gui.app %*
exit /b %ERRORLEVEL%