# Cómo poner el icono de **leona** en techamv.com (quitar la W de WordPress)

La **W** que ves suele ser:

1. **Favicon** del navegador (pestaña) → WordPress por defecto  
2. **Icono del sitio** (Site Icon) en el personalizador  
3. A veces el logo de la barra de admin de WP (solo si estás logueado)

## Archivos ya generados (carpeta `techamv-web/`)

| Archivo | Uso |
|---------|-----|
| `favicon.ico` | Clásico (pestaña del navegador) |
| `favicon-32.png` | 32×32 |
| `favicon-192.png` | Android / PWA |
| `favicon-180.png` | Apple touch |
| `leona-icon-512.png` | **Subir a WordPress** (Site Icon; WP recorta tamaños) |
| `assets/leona-icon-512.png` | Copia en assets |

---

## Opción A — WordPress (recomendado si techamv.com es WP)

1. Entra al **escritorio** de WordPress: `https://techamv.com/wp-admin`
2. Menú **Apariencia → Personalizar** (o **Apariencia → Editor** según tema)
3. **Identidad del sitio** / **Site Identity**
4. **Icono del sitio** / **Site Icon**
5. **Seleccionar archivo** → sube **`leona-icon-512.png`** (mínimo 512×512)
6. Ajusta el recorte (cuadrado, cara de la leona centrada)
7. **Publicar**

Eso sustituye la W en:

- pestaña del navegador  
- muchos menús móviles / favoritos  

Si no ves el cambio: vacía caché del hosting (LiteSpeed, Cloudflare, plugin de caché) y prueba en ventana privada.

### Si tu tema no tiene “Icono del sitio”

Plugin: **“Site Icon”** o cualquier SEO (Yoast / Rank Math suelen tener favicon),  
o sube `favicon.ico` a la **raíz del dominio** por FTP/File Manager:

```
public_html/favicon.ico
public_html/favicon-32.png
public_html/favicon-192.png
public_html/favicon-180.png
```

---

## Opción B — HTML estático (sin WP o página aparte)

En el `<head>` (ya añadido en `index.html` del repo):

```html
<link rel="icon" href="/favicon.ico" sizes="any">
<link rel="icon" type="image/png" href="/favicon-32.png" sizes="32x32">
<link rel="icon" type="image/png" href="/favicon-192.png" sizes="192x192">
<link rel="apple-touch-icon" href="/favicon-180.png">
```

Sube los archivos a la **raíz del sitio** (misma carpeta que `index.html`).

---

## Opción C — Solo quitar la W del admin (barra negra)

Eso es el logo de WordPress en el admin, no el favicon público.  
No hace falta para visitantes. Si quieres ocultarlo, es CSS de admin o un plugin de branding; lo importante para el público es el **Site Icon**.

---

## Checklist de subida (FTP / cPanel)

1. Conecta al hosting de **techamv.com**  
2. Ve a `public_html` (o la carpeta del dominio)  
3. Sube: `favicon.ico`, `favicon-32.png`, `favicon-192.png`, `favicon-180.png`, `leona-icon-512.png`  
4. En WP: Personalizar → Icono del sitio → `leona-icon-512.png`  
5. Ctrl+F5 o ventana privada para ver la leona en la pestaña  

---

*Icono generado para techamv · leona estilizada · 2026-07-23*
