# Commits — regla A

**Un commit = una carpeta de primer nivel.**  
El mensaje **nombra la carpeta** y **describe su contenido** (o el cambio).

## Formato

```
tipo(carpeta): resumen en una línea

Cuerpo:
- Qué es / qué contiene la carpeta
- Qué archivos toca el cambio
```

**Tipos:** `feat` · `fix` · `docs` · `chore` · `refactor`  
**Raíz del repo:** ámbito `raiz` (README, MAPA, FOLDERS, COMMITS, wiki si se considera bloque wiki → `wiki`).

## Ejemplos

```
docs(33x1): trato 1=repo civil por 33 años de paz — índice de la carpeta
feat(web-aleatorovix): motor JS sin Math.random + UI del ramo
docs(wiki): Home y ejes del monorepo
chore(techamv-aleatorovix-DEPLOY): pack FTP index + js + iconos
```

## Qué no hacer

- Un solo commit que mezcle `antipc/` + `Quijote/` + `33x1/` sin necesidad.  
- Mensajes tipo `Add files via upload` o `update stuff`.

Documento completo: [`COMMITS.md`](../COMMITS.md).
