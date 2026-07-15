@echo off
chcp 65001 >nul
cd /d "%~dp0"
echo === JUST RUN — Gemelos VMA ===
pip install -q numpy matplotlib 2>nul
cd gemelos
echo.
echo [1/3] gemelo_v6_fisic.py
python gemelo_v6_fisic.py
echo.
echo [2/3] hurto gravitatorio
cd ..\hurto-gravitatorio
python gemell_quijote_paper_rules.py
echo.
echo [3/3] zypyzape_twin_v4_8 (si existe)
cd ..\gemelos
if exist zypyzape_twin_v4_8_quijote.py (
  python zypyzape_twin_v4_8_quijote.py
) else (
  echo zypyzape_twin_v4_8_quijote.py no copiado — usar gemelo_v942.py
  python gemelo_v942.py
)
echo.
echo Listo. Revisa PNG y hurto-gravitatorio\quijote_results.csv
pause