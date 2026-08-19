# Aleatorovix 33×1 · Web del ramo ❤️🌹

Sitio público del **motor Aleatorovix** en el navegador (sin `Math.random`).

Autor: **Víctor Manzanares Alberola** · espiradesombra · VMA / 33×1  
Dominio de destino: [techamv.com](https://techamv.com)

## Estructura mínima para hosting (techamv)

```
index.html                 ← home (ramo, modo Princesa, Leonor)
aleatorovix/
  aleatorovix.js           ← motor (entropía → máscara Lila → bit flor)
```

Opcional: carpeta completa `aleatorovix/` (Python, C, demos, corpus).

## Probar en local

Abre `index.html` en el navegador (o sirve la carpeta con cualquier HTTP estático).

```bash
# Python
python -m http.server 8080
# luego http://localhost:8080
```

## API JS

```js
Aleatorovix.bit()        // 0 | 1
Aleatorovix.random()     // [0, 1)
Aleatorovix.randomInt(n) // [0, n)
Aleatorovix.ciclo()      // { decision, medida, k, bit, u, m }
```

## Python (paquete técnico)

```bash
cd aleatorovix
python aleatorovix.py
RUN_ALEATOROVIX.bat
```

## Subir a techamv.com (Nominalia / FTP)

1. Entra al panel de hosting (FTP o administrador de archivos).
2. Sube **`index.html`** a la **raíz pública** del dominio (suele ser `public_html`, `www` o `httpdocs`).
3. Crea la carpeta **`aleatorovix/`** en esa misma raíz.
4. Sube **`aleatorovix/aleatorovix.js`** dentro de esa carpeta.
5. Comprueba: `https://techamv.com/` y en la consola del navegador (F12) que no haya error 404 de `aleatorovix.js`.

**Importante:** si solo subes el HTML y no el `.js`, el ramo no usará Aleatorovix (fallará la carga del script).

Copia lista para subir: carpeta del Escritorio `techamv-aleatorovix-DEPLOY/`.

## Licencia / uso

Proyecto 33×1 · uso civil / educativo. Ver corpus en `aleatorovix/ALEATOROVIX.txt` y `LEEME.txt`.
