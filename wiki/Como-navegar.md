# Cómo navegar el monorepo

## Tres capas de documentación

| Capa | Dónde | Para qué |
|------|--------|----------|
| **Wiki** | `wiki/` (esta) | Visión humana, ejes, glosario |
| **Mapa / índice** | `MAPA.md`, `FOLDERS.md` | Lista de carpetas por tema o A–Z |
| **README de carpeta** | `carpeta/README.md` | Qué hay **dentro** de esa carpeta |

## Flujo recomendado

1. Lee [33x1](33x1) (el trato).  
2. Abre [Mapa-del-monorepo](Mapa-del-monorepo) y elige un eje.  
3. Entra en la carpeta del repo y lee su `README.md`.  
4. Ejecuta solo lo que el README / LEEME indiquen.

## Copia de trabajo canónica

| Ruta en el PC | Rol |
|---------------|-----|
| `Desktop\repo\claude-github\` | **Canónica** — monorepo completo, push aquí |
| `Desktop\web per llançar mobles\` | Borrador web → sincronizar a `web-aleatorovix/` |
| Carpetas sueltas del Escritorio | A menudo **más viejas** que el monorepo |

## Commits

Regla A: **un commit = una carpeta**, mensaje que describe contenido.  
Detalle: [Commits-regla-A](Commits-regla-A) y [`COMMITS.md`](../COMMITS.md).

## Regenerar índices

```bat
cd repo\claude-github
python tools\describe_folders.py
```

Eso actualiza `MAPA.md`, `FOLDERS.md` y READMEs auto; los README “de autor” se respetan salvo `--force`.
