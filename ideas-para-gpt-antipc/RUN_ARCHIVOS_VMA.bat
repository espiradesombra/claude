@echo off
chcp 65001 >nul
cd /d "%~dp0archivos-vma\codigo"
echo === ARCHIVOS-VMA — just run ===
echo.
echo [1] L(n) y m(n) — anexo E
python anexoE_L_m_script.py
echo.
echo [2] Calibrador MRAUV — anexo F
python anexoF_mrauv_calibrador.py n0=100000 dn=632
echo.
echo [3] Criba OpenMP (requiere gcc + OpenMP en WSL/MinGW):
echo     gcc -O3 -fopenmp anexoF_criba6kpm1_openmp.c -o criba ^&^& ./criba 1000000
echo.
echo Corpus: archivos-vma\CONCEPTO.txt
echo   cribas_cotas_vma.txt  — empezar aqui
pause