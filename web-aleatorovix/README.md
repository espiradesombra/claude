# web-aleatorovix (techamv / ramo)

Fuente canónica en Desktop:
`web per llançar mobles/aleatorovix/`

## Icono leona
- `favicon.ico`, `favicon-32.png`, … y `leona-icon-512.png`
- En el `<head>` de `index.html` (no quita easter eggs)

## Easter eggs (intactos)
- Escribe **leonor** (teclado) → panel científico
- Escribe **princesa de asturias su majestad de españa** → modo Princesa (3301, mód, doble mód)
- Escribe **normal** → sale del modo Princesa

## Entropía / “red”
En `aleatorovix.js`, `red_x` mezcla:
1. jitter de CPU local (como siempre)
2. **opcional** `navigator.connection` (rtt/downlink) + online/offline si el navegador lo da

No depende de red: si no hay API de red, sigue funcionando.

## Subir a techamv.com
Sube carpeta deploy `techamv-aleatorovix-DEPLOY/` o esta carpeta con `index.html` + `aleatorovix.js` + favicons en la raíz del site.
