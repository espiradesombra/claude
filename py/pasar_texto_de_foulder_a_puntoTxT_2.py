"""
pasar_texto_de_foulder_a_puntoTxT.py
Autor: VMA
Convierte todos los .docx y .pdf de una carpeta (y subcarpetas) a .txt
Uso:
    python pasar_texto_de_foulder_a_puntoTxT.py
    python pasar_texto_de_foulder_a_puntoTxT.py /ruta/a/carpeta
"""

import sys
import os
import glob

# ── dependencias ──────────────────────────────────────────────────────────────
try:
    from docx import Document
except ImportError:
    print("Instalando python-docx...")
    os.system("pip install python-docx --break-system-packages -q")
    from docx import Document

try:
    import fitz  # pymupdf
except ImportError:
    print("Instalando pymupdf...")
    os.system("pip install pymupdf --break-system-packages -q")
    import fitz

# ── carpeta de trabajo ────────────────────────────────────────────────────────
carpeta = sys.argv[1] if len(sys.argv) > 1 else "."
carpeta = os.path.abspath(carpeta)
print(f"\n📂 Carpeta: {carpeta}\n")

convertidos = 0
errores = 0

# ── recorrer archivos ─────────────────────────────────────────────────────────
for path in sorted(glob.glob(os.path.join(carpeta, "**", "*"), recursive=True)):

    ext = os.path.splitext(path)[1].lower()
    if ext not in (".docx", ".pdf"):
        continue

    out = os.path.splitext(path)[0] + ".txt"
    nombre = os.path.relpath(path, carpeta)

    try:
        if ext == ".docx":
            doc = Document(path)
            texto = "\n".join(p.text for p in doc.paragraphs)

        elif ext == ".pdf":
            doc = fitz.open(path)
            texto = "\n".join(page.get_text() for page in doc)
            doc.close()

        if not texto.strip():
            print(f"⚠️  vacío  →  {nombre}")
            continue

        with open(out, "w", encoding="utf-8") as f:
            f.write(texto)

        kb = os.path.getsize(out) // 1024
        print(f"✓  {nombre}  →  .txt  ({kb} KB)")
        convertidos += 1

    except Exception as e:
        print(f"✗  {nombre}  →  ERROR: {e}")
        errores += 1

# ── resumen ───────────────────────────────────────────────────────────────────
print(f"\n{'─'*50}")
print(f"✅ Convertidos: {convertidos}   ❌ Errores: {errores}")
print(f"{'─'*50}\n")
