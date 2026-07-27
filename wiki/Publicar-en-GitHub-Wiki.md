# Publicar en la pestaña Wiki de GitHub

El monorepo ya tiene la wiki en la carpeta **`wiki/`** (versionada con el código).  
La pestaña **Wiki** de GitHub (`…/claude/wiki`) es un **repositorio git aparte** y solo aparece con contenido cuando existe la primera página.

## Estado actual

- Repo: `has_wiki: true`  
- Clone de `claude.wiki.git`: vacío / “not found” hasta crear la 1ª página en la web o con credenciales.

## Opción A — Usar solo `wiki/` del monorepo (recomendada ya)

1. Lee las páginas en:  
   https://github.com/espiradesombra/claude/tree/main/wiki  
2. Enlaces desde el README raíz → sección Wiki.  
3. Commits: `docs(wiki): …` (regla A, carpeta `wiki`).

## Opción B — Activar la pestaña Wiki de GitHub

1. En GitHub: **Settings → Features → Wikis** (ya suele estar ON).  
2. Abre **Wiki** → *Create the first page* (puede ser un “Home” mínimo).  
3. En local:

```bat
git clone https://github.com/espiradesombra/claude.wiki.git
cd claude.wiki
copy /Y ..\claude-github\wiki\*.md .
git add *.md
git commit -m "docs(wiki): import monorepo wiki Home and axes"
git push origin master
```

(En algunos repos la rama del wiki es `master`, no `main`.)

4. Sidebar: el archivo `_Sidebar.md` lo usa GitHub Wiki automáticamente.

## Opción C — Pedir al agente

Cuando la primera página exista (o haya token con permiso), el agente puede clonar `claude.wiki.git`, copiar `wiki/*.md` y hacer push a la pestaña Wiki.

## Nombres de página

GitHub Wiki: `Mapa-del-monorepo.md` → URL `…/wiki/Mapa-del-monorepo`.  
Mantener los mismos nombres de archivo que en `wiki/` del monorepo.
