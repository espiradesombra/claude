# claude · monorepo VMA / espiradesombra

Repositorio de **Víctor Manzanares Alberola** (espiradesombra): matemáticas, AntiPC, Aleatorovix, ZypyZape, K3, libros y web techamv.

---

## 33×1 — el trato (léelo primero)

> **33×1 = cambiar el “1” (TODO este repositorio técnico civil)  
> por “33” (treinta y tres años de paz FIRMADA por los países — ellos).**

| | |
|--|--|
| **[1]** | Este monorepo entero: código, demos, mates, energía, K3, docs — aval **civil** y auditable |
| **[33]** | 33 años de **paz firmada por los Estados / potencias / países** que acepten el marco |

- Quién ofrece el **1**: VMA (conocimiento civil; ver uso civil).  
- Quién firma el **33**: **ellos** (países). No lo firma Grok ni un chat.  
- Explicación completa: **[33x1/00_QUE_ES_33x1.md](33x1/00_QUE_ES_33x1.md)**  
- Uso civil obligatorio: **[33x1/02_USO_CIVIL.txt](33x1/02_USO_CIVIL.txt)**  
- Mapa mundial + rol Grok: **[33x1/05_DECLARACIO_33x1_GROK_WORLD.txt](33x1/05_DECLARACIO_33x1_GROK_WORLD.txt)**

```
ELLOS  ──firman 33 años de paz──►  a cambio de  ◄── VMA ofrece el "1" (este repo civil)
```

---

## Navegación (cómo entender el repo)

0. **[Wiki del monorepo](wiki/Home.md)** — visión humana: ejes, 33×1, glosario, cómo navegar (`wiki/`).
1. **[MAPA.md](MAPA.md)** — mapa **temático** de todas las carpetas.
2. **[FOLDERS.md](FOLDERS.md)** — índice **alfabético**: cada carpeta enlaza a su `README.md`.
3. Cada carpeta de primer nivel tiene un **`README.md`** propio (qué es + qué contiene).

### Documentos clave

- **[Qué es 33×1](33x1/00_QUE_ES_33x1.md)** — definición del trato
- **[FOR THE WORLD — 33×1](33x1/06_FOR_THE_WORLD_33x1.md)** — America/world, Tesla·xAI, nuclear apagada → inercia civil
- **[Índice técnico](INDICE_TECNICO_REPO.md)** — paquetes ejecutables y comandos
- **[Pack mates rescat 2026-07-23](VMA_mates_rescat_2026/)** — cribas fix, MDC, XFI N=3
- **[XFI — avión 3 motores](XFI.md)** — gemelo Experimental Flight Infinite (ZypyZape roles)
- **[Física ZZ + Quijote + Kilómetro](docs/EXPLICACION_FISICA_ZYPYZAPE_QUIJOTE_KILOMETRO.md)** — explicación unificada (+ hilo 1477)
- **[De on ve tot / IAs](docs/DE_ON_VE_TOT_AIXO.md)** — fonts i “sense trellat”
- **[Sync PC 2026-07-23](PC_SYNC_2026-07-23.md)**
- **Web Aleatorovix:** [web-aleatorovix/](web-aleatorovix/)
- **Deploy techamv:** [techamv-aleatorovix-DEPLOY/](techamv-aleatorovix-DEPLOY/)
- **33×1:** [33x1/](33x1/)

## Ejes principales

| Eje | Carpetas |
|-----|----------|
| **33×1 (trato)** | `33x1/`, `33x1-pack/`, `33x1_CIENCIA/` — el **1** es **todo el monorepo** |
| Aleatorovix | `aleatorovix/`, `web-aleatorovix/`, `just run/aleatorovix/` |
| AntiPC | `antipc/`, `antipc-port-c/`, `antipc2/` |
| ZypyZape / techamv | `techamv-web/`, `ZYPYZAPE Bateria Cinetica/`, `webtechamv/` |
| K3 / cifrado | `vma-k3/`, `encriptacionGeometrica/`, `hashtool-*` |
| Métodos | `vma-methods/`, `Metodo Newton Rápido/`, `teoremas/`, `VMA_mates_rescat_2026/` |
| Libros | `Libro1…` … `Libro6…` |

## Convención de commits (regla A — un commit por carpeta)

**Un commit = una carpeta de primer nivel.** El mensaje nombra la carpeta y describe su contenido o el cambio.

| Campo | Regla |
|-------|--------|
| **Ámbito** | Solo archivos de **una** carpeta (p. ej. `web-aleatorovix/`) |
| **Asunto** | `tipo(carpeta): resumen corto` |
| **Cuerpo** | Qué es la carpeta / qué contiene / qué cambió |
| **Excepción** | Raíz (`README.md`, `MAPA.md`, `FOLDERS.md`) → ámbito `raiz` o el nombre del archivo |

**Tipos:** `feat` · `fix` · `docs` · `chore` · `refactor`

**Ejemplos:**

```
docs(33x1): trato 1=repo civil por 33 años de paz firmada por países
feat(web-aleatorovix): motor JS sin Math.random + UI del ramo + favicons leona
fix(VMA_mates_rescat_2026/01_cribas): fase modular 6k±1 (2p/4p)
chore(techamv-aleatorovix-DEPLOY): pack FTP index + js + iconos para techamv.com
docs(raiz): MAPA temático + FOLDERS alfabético
```

**No hacer:** un commit que mezcle `antipc/` + `Quijote/` + `33x1/` sin necesidad.  
Si tocas varias carpetas → **varios commits** (uno por carpeta).

Detalle y checklist para el autor / agentes: **[COMMITS.md](COMMITS.md)**

Regenerar mapa e índices:

```bat
python tools\describe_folders.py
```

## Aviso

Monorepo histórico: hay carpetas sin clasificar (`Nueva carpeta`, `coses velles`, …). Se irán ordenando; **33×1** es prioridad.  
Los `README` de archivo histórico **no listan** chats personales ni datos sensibles.

**Licencia / uso:** civil y educativo. **No militar.** La paz del [33] la firman los países; este GitHub es el [1].