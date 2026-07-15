@echo off
chcp 65001 >nul
cd /d "%~dp0"
echo Compilando metodos_simples_vma.tex ...
pdflatex -interaction=nonstopmode metodos_simples_vma.tex
pdflatex -interaction=nonstopmode metodos_simples_vma.tex
if exist metodos_simples_vma.pdf (
  echo.
  echo OK: metodos_simples_vma.pdf
) else (
  echo.
  echo ERROR: instala MiKTeX o TeX Live y pdflatex en PATH
)
pause