# Aleatorovix y web (techamv)

Motor de entropía / decisiones **sin `Math.random`**: máscara Lila, criba, MDC. UI del ramo ❤️🌹 para [techamv.com](https://techamv.com).

## Carpetas

| Carpeta | Contenido |
|---------|-----------|
| `aleatorovix/` | Organismo Python/C, demos, benchmarks, corpus |
| `web-aleatorovix/` | Sitio completo: HTML, JS, favicons leona, paquete técnico |
| `techamv-aleatorovix-DEPLOY/` | **Pack FTP mínimo**: index + js + iconos + guía COMO_SUBIR |
| `techamv-web/` | Sitio histórico ZypyZape/TechAMV (HTML/PHP) |
| `webtechamv/` | Hubs K3 / licencia (React TSX) |
| `desktop-snapshot/` | Copias sueltas del Escritorio (no canónicas) |

## Publicar en techamv.com

1. Usar solo `techamv-aleatorovix-DEPLOY/`.  
2. Subir a `public_html` (o equivalente) por FTP.  
3. Debe existir `aleatorovix/aleatorovix.js` junto al `index.html`.  
4. Comprobar en F12 que no hay 404 del script.

Guía: `techamv-aleatorovix-DEPLOY/COMO_SUBIR_TECHAMV.txt`.

## API JS (resumen)

```js
Aleatorovix.bit()
Aleatorovix.random()
Aleatorovix.randomInt(n)
Aleatorovix.ciclo()
```
