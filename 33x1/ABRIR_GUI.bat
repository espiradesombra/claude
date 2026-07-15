@echo off
chcp 65001 >nul
cd /d "C:\Users\cuent\Desktop\repo\vma-methods"
if exist vma-gui.cmd (
    call vma-gui.cmd
) else (
    echo No encontrado repo\vma-methods\vma-gui.cmd
    pause
)