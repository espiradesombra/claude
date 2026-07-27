# Convención de commits — monorepo `claude`

**Regla A (activa):** **un commit por carpeta** de primer nivel.  
El commit debe **nombrar la carpeta** y **describir su contenido** (o el cambio concreto).

Repo: [espiradesombra/claude](https://github.com/espiradesombra/claude)  
Autor: Víctor Manzanares Alberola · espiradesombra · VMA / 33×1

---

## Formato

```
tipo(carpeta): resumen en una línea

Cuerpo (opcional pero recomendado):
- Qué es esta carpeta
- Qué archivos / subcarpetas toca el cambio
- Por qué (si no es obvio)
```

### Tipos

| Tipo | Uso |
|------|-----|
| `feat` | Función o material nuevo |
| `fix` | Corrección |
| `docs` | README, mapas, textos, papeles |
| `chore` | Empaquetado, zip, deploy, limpieza |
| `refactor` | Reordenar sin cambiar comportamiento |

### Ámbito `(carpeta)`

- Nombre de la carpeta de **primer nivel**: `web-aleatorovix`, `33x1`, `antipc`, …
- Subcarpeta importante: `VMA_mates_rescat_2026/01_cribas`
- Archivos de la **raíz** del repo: `raiz` (o `MAPA`, `FOLDERS`, `README`)

### Ejemplos buenos

```
docs(33x1): uso civil + definición del trato 33×1
feat(aleatorovix): organismo entropía + máscara Lila + demos
chore(techamv-aleatorovix-DEPLOY): pack mínimo FTP para techamv.com
fix(antipc-port-c): criba híbrida port C + inventarios v04–v14
docs(raiz): MAPA temático y convención un-commit-por-carpeta
```

### Ejemplos malos

```
update stuff
Add files via upload
varios cambios
fix: todo el monorepo
```

---

## Flujo de trabajo

1. Trabajar **solo** en una carpeta (o hacer stash / commits parciales).
2. Actualizar el `README.md` de esa carpeta si el contenido cambia de significado.
3. Si hay carpetas nuevas o renombres:  
   `python tools\describe_folders.py`  
   y commit aparte de `docs(raiz)` para MAPA/FOLDERS si cambian.
4. `git add` **solo** rutas de esa carpeta.
5. Commit con el formato de arriba.
6. Si hay otra carpeta → repetir (otro commit).
7. `git push origin main` cuando toque publicar.

### Staging parcial (varias carpetas tocadas)

```bat
git add web-aleatorovix/
git commit -m "feat(web-aleatorovix): ..."

git add techamv-aleatorovix-DEPLOY/
git commit -m "chore(techamv-aleatorovix-DEPLOY): ..."
```

---

## Qué necesita el agente (Grok / humano coautor)

Para aplicar la regla A sin adivinar:

| Necesito de ti | Ejemplo |
|----------------|---------|
| **Carpeta** (una) | `web-aleatorovix` o `33x1` |
| **Qué hacer** | “sube el favicon”, “arregla criba”, “documenta X” |
| **Origen en el PC** si no está en el monorepo | `Desktop\web per llançar mobles\…` |
| **¿Push a GitHub?** | “solo commit” / “sube” |
| **Prioridad 33×1** si aplica | “esto es del 1 civil” / “no tocar 33x1” |

### No hace falta que digas

- El hash del commit
- Cada archivo suelto si la carpeta está clara
- “haz un commit gigante de todo” (lo partiremos por carpeta)

### No subir

- Secretos (FTP, tokens, `.env`)
- Chats personales / teléfonos / NIF
- Vídeos enormes / cracks / basura del Escritorio

---

## Copia de trabajo canónica

| Ruta local | Rol |
|------------|-----|
| `Desktop\repo\claude-github\` | **Canónica** — monorepo completo, push aquí |
| `Desktop\web per llançar mobles\` | Borrador web; sincronizar hacia `web-aleatorovix/` antes de commit |

---

*Regla acordada en sesión 2026-07-23: opción A — un commit por carpeta.*
