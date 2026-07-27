# techamv-aleatorovix-DEPLOY

**Pack mínimo para publicar el ramo Aleatorovix en [techamv.com](https://techamv.com)** (FTP / panel del hosting).

Autor: Víctor Manzanares Alberola · espiradesombra · VMA / 33×1

## Qué contiene esta carpeta

| Ruta | Contenido |
|------|-----------|
| `index.html` | Página del ramo ❤️🌹 (UI, modos, número del destino) |
| `aleatorovix/aleatorovix.js` | Motor en el navegador **sin `Math.random`** (máscara Lila + entropía) |
| `favicon-*.png` / `favicon.ico` | Iconos del sitio (leona) 16–512 px |
| `leona-icon-512.png` | Icono grande / WordPress site icon |
| `COMO_SUBIR_TECHAMV.txt` | Pasos Nominalia/FTP para subir este pack |
| `README.md` | Este archivo |

## Cómo se usa

1. Copia de seguridad del `index.html` actual en el hosting.
2. Sube **todo el contenido** de esta carpeta a la raíz pública (`public_html` / `www` / `httpdocs`).
3. Debe quedar: `…/index.html` y `…/aleatorovix/aleatorovix.js`.
4. Abre https://techamv.com y comprueba en F12 que no hay 404 del `.js`.

## Relación con el monorepo

- Código completo y demos: [`../web-aleatorovix/`](../web-aleatorovix/)
- Organismo Python/C: [`../aleatorovix/`](../aleatorovix/)
- Mapa del monorepo: [`../MAPA.md`](../MAPA.md)

**No** subas secretos FTP al Git. Este pack es solo estáticos públicos.
