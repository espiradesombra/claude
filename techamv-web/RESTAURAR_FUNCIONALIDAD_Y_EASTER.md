# Restaurar funcionalidad + easter eggs + icono leona

## Qué pasó
Al añadir el favicon se usó/publicó un `index.html` **más corto** (sin todo el contacto/analytics)
o se subió solo iconos y se tocó el site en WP. Los easter eggs viven sobre todo en la **404 Data Room**.

## Qué hay ahora (arreglado en el repo)

| Archivo | Contenido |
|---------|-----------|
| `index.html` | Home **completa** (chat, contacto PHP, analytics, NDA/Data Room) + favicon leona |
| `404.html` | Data Room con easter eggs (**4 clics** + **F11**) + favicon leona |
| `favicon.*` / `leona-icon-512.png` | Iconos (no quitan JS) |

Backup del index delgado: `index-thin-solo-favicon-backup.html`

## Cómo subir a techamv.com (cPanel / FTP)

1. Sube a la **raíz del sitio** (o la carpeta que sirva el HTML estático de ZypyZape):
   - `index.html` (el nuevo completo)
   - `404.html` (importante para easter eggs)
   - `favicon.ico`, `favicon-32.png`, `favicon-192.png`, `favicon-180.png`
   - `leona-icon-512.png` (para WP Site Icon)
   - `chat.php`, `enviar_correo.php` si no están ya en `/api/` o donde apunten

2. En WordPress, si la home es una **página estática** o un tema:
   - No sustituyas toda la instalación WP por un HTML a ciegas.
   - Si la landing es un HTML custom en el dominio:
     - Home = `index.html`
     - Error 404 del servidor / página 404 de WP debe apuntar al Data Room HTML
   - En WP: **Ajustes → Lectura** o redirección 404 del hosting a `/404.html`
   - En WP: **Apariencia → Personalizar → Icono del sitio** = `leona-icon-512.png`

3. Prueba:
   - Home: chat y contacto funcionan
   - `/404` o URL inventada: Data Room + 4 clics en logo + escribe F11 en factor
   - Pestaña del navegador: leona, no W

## Easter eggs (recordatorio)

- **404 / Data Room**: 4 clics en el gatillo secreto; input **F11** desclasifica memorándum
- Enlaces home → `https://techamv.com/equipo/404` (Data Room)

Si tras subir aún falta algo, di qué URL exacta falla.
