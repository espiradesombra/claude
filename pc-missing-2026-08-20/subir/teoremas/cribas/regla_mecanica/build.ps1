# build.ps1 — ReglaMecanicaUniversal.dll + test
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $here

function Find-Gpp {
    @("g++", "x86_64-w64-mingw32-g++", "clang++") | ForEach-Object {
        if (Get-Command $_ -ErrorAction SilentlyContinue) { return $_ }
    }
    $paths = @(
        "C:\msys64\mingw64\bin\g++.exe",
        "C:\MinGW\bin\g++.exe",
        "C:\Program Files\LLVM\bin\clang++.exe"
    )
    foreach ($p in $paths) {
        if (Test-Path $p) { return $p }
    }
    return $null
}

$cpp = Find-Gpp
if (-not $cpp) {
    Write-Host "No se encontró g++/clang++. Usa regla_mecanica.py para probar la lógica."
    python regla_mecanica.py
    exit 1
}

Write-Host "Compilador: $cpp"
& $cpp -O2 -std=c++17 test_regla_mecanica.cpp -o test_regla_mecanica.exe -lm
if ($LASTEXITCODE -eq 0) {
    & ".\test_regla_mecanica.exe"
}

& $cpp -shared -O2 -std=c++17 -DBUILDING_DLL dll_export.cpp -o ReglaMecanicaUniversal.dll -lm
if ($LASTEXITCODE -eq 0) {
    Write-Host "OK: ReglaMecanicaUniversal.dll"
} else {
    Write-Host "Fallo DLL (revisar exports)"
}