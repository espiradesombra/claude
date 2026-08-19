# Hash K3 — portables y elimina-duplicados

## ¿Hay portables?

Sí, **tres herramientas** en el pack C (`k3hash`):

| Tool | Qué hace |
|------|----------|
| **`k3dedup`** | Duplicados **exactos** (mismo tamaño + mismo hash K3) |
| **`k3simil`** | Similitud / plagio parcial (shingles + Jaccard) |
| **`k3search`** | Índice + búsqueda de contenido en textos |

Código fuente (copia canónica en el monorepo):

`hashtool-work/k3hash/k3hash/`  
también en `inbox/subir/VMA-Hash-Pack-v1-staging/k3hash/`

**Hoy en este PC:** no hay `k3dedup.exe` precompilado (solo `k3hash.dll` en AntiPC y el `.c`).  
Para usar **ya**, sin CMake:

```powershell
cd C:\Users\cuent\Desktop\repo\claude-github\hashtool-work
python k3dedup.py C:\ruta\a\escanear
```

Ejemplo de prueba (a.txt ≈ b.txt, c.txt distinto):

```powershell
python k3dedup.py ..\docs\integracion-antipc\k3dedup_demo
```

Opciones útiles:

```powershell
python k3dedup.py C:\Datos --ext .pdf,.docx,.py
python k3dedup.py C:\Datos --delete-report   # sugiere qué borrar (no borra)
# por tubería:
Get-ChildItem -Recurse -File C:\Datos | ForEach-Object FullName | python k3dedup.py -
```

No borra ficheros: solo imprime grupos. Tú decides.

---

## Compilar los 3 portables C (cuando tengas CMake + compilador)

```powershell
cd hashtool-work\k3hash\k3hash
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

Salida típica: `build\Release\k3dedup.exe`, `k3simil.exe`, `k3search.exe`, `k3hash.dll`.

Uso Windows del `.exe`:

```bat
dir /s /b C:\mi\carpeta | build\Release\k3dedup.exe
```

Linux/macOS:

```bash
find /mi/carpeta -type f | ./k3dedup
find /mi/carpeta -name "*.txt" | ./k3simil 0.30
```

Detalle completo: `hashtool-work/k3hash/k3hash/README.md`

---

## AntiPC (plugin, no portable suelto)

`antipc/plugins/k3_dedup_plugin.py` — mismo criterio (tamaño + digest) dentro del runtime de plugins. No es un `.exe` de escritorio.

---

## GUI Windows (recomendado): EliminaDuplicadosK3

Ejecutable listo (hash K3 + borrado con confirmación `BORRAR`):

- `hashtool-work/bin/EliminaDuplicadosK3.exe`
- Fuente: `hashtool-work/elimina_duplicados_k3.py`
- CLI solo-listar: `python hashtool-work/k3dedup.py C:\carpeta`

Doble confirmación; nunca borra todas las copias de un grupo.

