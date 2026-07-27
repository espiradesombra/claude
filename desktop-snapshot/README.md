# desktop-snapshot

**Instantánea** de archivos sueltos del Escritorio que conviene tener en el monorepo sin mezclarlos con el pack web canónico.

Autor: VMA · referencia de trabajo, no siempre la versión “de producción”.

## Qué contiene esta carpeta

| Archivo | Contenido |
|---------|-----------|
| `index-aleatorovix-techamv.html` | Copia HTML del ramo (referencia) |
| `aleatorovix.js` | Copia del motor (puede ir **detrás** de `web-aleatorovix/`) |
| `COMO_SUBIR_TECHAMV.txt` | Guía FTP / Nominalia para techamv.com |
| `MAPA_VIEJO_Y_V.txt` | Mapa viejo del PC / notas de organización |
| `Nuevo Documento de texto (32).txt` | Notas de sesión / borrador |

## Canónico (usar esto para publicar)

| Uso | Carpeta |
|-----|---------|
| Web completa + demos | [`../web-aleatorovix/`](../web-aleatorovix/) |
| Pack FTP mínimo | [`../techamv-aleatorovix-DEPLOY/`](../techamv-aleatorovix-DEPLOY/) |
| Mapa monorepo | [`../MAPA.md`](../MAPA.md) |

Si `aleatorovix.js` aquí y en `web-aleatorovix/` no coinciden, **gana** `web-aleatorovix/`.
