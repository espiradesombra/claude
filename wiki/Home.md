# Wiki del monorepo `claude`

**Autor:** Víctor Manzanares Alberola (espiradesombra) · VMA  
**Repo:** [github.com/espiradesombra/claude](https://github.com/espiradesombra/claude)  
**Dominio web (ramo):** [techamv.com](https://techamv.com)

Esta wiki **describe el repositorio**: qué es, por dónde empezar, qué hay en cada eje y cómo no perderse en ~120 carpetas de primer nivel.

---

## En una frase

Monorepo civil de matemáticas (cribas, MDC, Newton, Goldbach…), AntiPC, Aleatorovix, cifrado K3, energía (ZypyZape / Quijote / Kilómetre), libros VMA y el trato **33×1**.

---

## 33×1 (léelo antes que el código)

> **33×1 = cambiar el “1” (TODO este monorepo técnico civil)  
> por “33” (33 años de paz firmada por los países).**

| Pieza | Quién | Qué |
|-------|--------|-----|
| **[1]** | VMA ofrece | Todo el repo civil y auditable |
| **[33]** | Los países firman | 33 años de paz activa |

- Definición: [`33x1/00_QUE_ES_33x1.md`](../33x1/00_QUE_ES_33x1.md)
- Uso civil: [`33x1/02_USO_CIVIL.txt`](../33x1/02_USO_CIVIL.txt)
- Mundo: [`33x1/06_FOR_THE_WORLD_33x1.md`](../33x1/06_FOR_THE_WORLD_33x1.md)
- Página wiki: [33x1](33x1)

**Uso:** civil y educativo. **No militar.**

---

## Qué usar YA (criterio del monorepo)

| Prioridad | Carpeta / página | Por qué |
|:---------:|------------------|---------|
| **1** | **[Demo oficial Sale-It](Demo-oficial-Sale-It)** → `sale-it/` | Producto con bats, K3 y benchmark A–E medido |
| **2** | `techamv-aleatorovix-DEPLOY/` | Pack FTP mínimo para techamv.com |
| **3** | `vma-methods/` | Mates por CLI en un minuto |
| **4** | `VMA_mates_rescat_2026/` | Pack mates limpio con tests |
| **5** | Resto del monorepo | Corpus / histórico — no es la demo |

**No** priorizar rescatar “monton viejo” ni duplicados de `subir/`: casi todo lo usable **ya está** en GitHub; el ROI está en **demostrar y empaquetar**, no en re-subir zips.

---

## Por dónde empezar (ruta de 5 minutos)

| Orden | Ir a | Para qué |
|------:|------|----------|
| 0 | [Demo oficial Sale-It](Demo-oficial-Sale-It) | Correr algo con números |
| 1 | [33x1](33x1) | Entender el trato |
| 2 | [Mapa-del-monorepo](Mapa-del-monorepo) | Ver ejes y carpetas |
| 3 | [Como-navegar](Como-navegar) | README por carpeta, FOLDERS, commits |
| 4 | Eje que te interese | [AntiPC](AntiPC) · [Aleatorovix-y-web](Aleatorovix-y-web) · [Matematicas](Matematicas) · [Energia-ZypyZape-Quijote](Energia-ZypyZape-Quijote) · [K3-y-cifrado](K3-y-cifrado) · [Libros-VMA](Libros-VMA) |
| 5 | [Commits-regla-A](Commits-regla-A) | Un commit = una carpeta |

En el repo (raíz):

- [`MAPA.md`](../MAPA.md) — mapa temático  
- [`FOLDERS.md`](../FOLDERS.md) — índice alfabético  
- [`COMMITS.md`](../COMMITS.md) — convención de commits  
- Cada carpeta → su `README.md`

---

## Ejes del monorepo

| Eje | Wiki | Carpetas típicas |
|-----|------|------------------|
| Trato 33×1 | [33x1](33x1) | `33x1/`, `33x1-pack/`, `33x1_CIENCIA/` |
| Mates / métodos | [Matematicas](Matematicas) | `vma-methods/`, `teoremas/`, `VMA_mates_rescat_2026/`, Goldbach, Newton… |
| AntiPC | [AntiPC](AntiPC) | `antipc/`, `antipc-port-c/`, `sale-it/` |
| Aleatorovix + web | [Aleatorovix-y-web](Aleatorovix-y-web) | `aleatorovix/`, `web-aleatorovix/`, `techamv-aleatorovix-DEPLOY/` |
| K3 / cifrado | [K3-y-cifrado](K3-y-cifrado) | `vma-k3/`, `encriptacionGeometrica/`, `hashtool-*` |
| Energía | [Energia-ZypyZape-Quijote](Energia-ZypyZape-Quijote) | `ZYPYZAPE…`, `Quijote/`, `kilometre…`, `docs/` |
| Libros | [Libros-VMA](Libros-VMA) | `Libro1` … `Libro6` |
| Empaquetado | [Paquetes-y-ejecucion](Paquetes-y-ejecucion) | `just run/`, `sale-it/`, `bin/` |

---

## Qué no es esta wiki

- No sustituye el código ni los papers.
- No es el hosting de techamv.com (eso es FTP con el pack DEPLOY).
- No lista chats personales ni secretos.

---

## Mantenimiento

- Contenido versionado en la carpeta [`wiki/`](../wiki/) del monorepo (commit `docs(wiki): …`).
- La pestaña **Wiki** de GitHub se puede sincronizar más adelante (ver [Publicar-en-GitHub-Wiki](Publicar-en-GitHub-Wiki)).
- Regenerar mapas de carpetas: `python tools/describe_folders.py`
