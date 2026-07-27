# web-aleatorovix

**Sitio web del ramo Aleatorovix** para el monorepo y para techamv.com: motor JS en el navegador **sin `Math.random`**, UI del ramo ❤️🌹, favicons leona y paquete técnico Python/C.

Autor: Víctor Manzanares Alberola · espiradesombra · VMA / 33×1

## Qué contiene esta carpeta

| Ruta | Contenido |
|------|-----------|
| `index.html` | Página principal del ramo (despliegue / local) |
| `aleatorovix.js` | Motor JS en la raíz (mismo núcleo que la subcarpeta) |
| `aleatorovix/` | Paquete: `aleatorovix.js`, `index.html`, Python, C, demos, corpus, **favicons** |
| `aleatorovix/nucleo/` | Módulos: entropía, máscara Lila, criba, MDC memoria, acciones |
| `aleatorovix/demos/` | Demos primos / MDC / organismo |
| `aleatorovix/benchmarks/` | Benches y resultados |
| `favicon-*` / `leona-icon-512.png` | Iconos del sitio |
| `aleatorovix_san_valent_n*.html`, `indexnv.html`, … | Variantes / historial HTML |
| `ALEATOROVIX.txt`, `README.md` | Corpus y esta guía |

## Probar en local

```bat
cd web-aleatorovix
python -m http.server 8080
```

Abre `http://localhost:8080` (o el `index.html` directamente).

## API JS (resumen)

```js
Aleatorovix.bit()        // 0 | 1
Aleatorovix.random()     // [0, 1)
Aleatorovix.randomInt(n) // [0, n)
Aleatorovix.ciclo()      // { decision, medida, k, bit, u, m }
```

## Deploy a techamv.com

Usar el pack mínimo: [`../techamv-aleatorovix-DEPLOY/`](../techamv-aleatorovix-DEPLOY/)  
(no hace falta subir todo el corpus Python/C al hosting).

## Enlaces

- Mapa monorepo: [`../MAPA.md`](../MAPA.md)
- 33×1: [`../33x1/`](../33x1/)
- Organismo canónico extra: [`../aleatorovix/`](../aleatorovix/)
