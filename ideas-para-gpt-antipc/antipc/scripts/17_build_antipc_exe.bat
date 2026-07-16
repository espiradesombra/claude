@echo off
chcp 65001 >nul
cd /d "%~dp0.."
echo === Build antipc.exe (PyInstaller) ===
echo.

python -m pip install pyinstaller --quiet
if errorlevel 1 (
    echo ERROR: pip install pyinstaller
    exit /b 1
)

if not exist "dist" mkdir dist
if not exist "build\work" mkdir build\work

python -m PyInstaller --noconfirm --clean --distpath dist --workpath build\work build\antipc.spec
if errorlevel 1 (
    echo ERROR: pyinstaller fallo
    exit /b 1
)

echo.
echo OK: dist\antipc.exe
dist\antipc.exe version
exit /b 0