---
name: conceptual-map-of-chat
description: >
  creativeTxt — responde escribiendo un único archivo .txt artístico (mapa
  conceptual del chat). ASCII art denso, cajas, ramas, pulso temporal. Usar
  cuando el usuario pida creativeTxt, arthacher, mapa conceptual, responder
  en txt, arte en texto, o ejecute /conceptual-map-of-chat. El chat solo
  confirma la ruta; el contenido va en el archivo.
metadata:
  short-description: "creativeTxt — mapa artístico en .txt"
---

# conceptual-map-of-chat (creativeTxt)

Transforma la conversación en un **único archivo `.txt` artístico**: mapa conceptual visual del chat.

## Regla de oro

**No respondas el contenido en el chat.** El chat recibe solo:

```
🎨 creativeTxt → <ruta/absoluta/al/archivo.txt>
Abre el archivo para ver el mapa.
```

Máximo 3 líneas en chat. Todo el arte va al `.txt`.

## Antes de escribir

1. Lee `references/arthacher.txt` — manual de estilo.
2. Repasa el hilo del chat: tema, bifurcaciones, ideas clave, tono.
3. Elige nombre: `conceptual-map-YYYYMMDD-HHMM.txt` o `creative-<tema>.txt`.
4. Carpeta por defecto: `C:\Users\cuent\Desktop\creative-maps\` (crear si no existe).
5. Si el usuario da otra ruta, usarla.

## Contenido obligatorio del .txt

```
╔══════════════════════════════════════════════════╗
║  TÍTULO DEL MAPA                    fecha/hora   ║
╚══════════════════════════════════════════════════╝

[1] NÚCLEO — idea central en caja grande

[2] RAMAS — 3–7 nodos con líneas ┌ │ └ ─ →

[3] PULSO — línea temporal del chat (inicio → ahora)

[4] DETALLE — citas o conceptos del usuario en mini-cajas

[5] SÍNTESIS — cierre visual (diamante, estrella, espiral ASCII)

[6] PIE — firma: ◆ conceptualMapOfChat ◆ + semilla visual única
```

## Estilo artístico

- **Alta densidad**: bordes dobles `╔═╗`, sombras `█▓▒░`, conectores `├─┤`.
- **Forma = significado**: árbol para jerarquía, espiral para evolución, puente `═══` para enlaces.
- **Unicode permitido**: → ↑ ↓ ➜ ★ ✦ ◇ ╭ ╮ (mantener legible en Notepad).
- **Sin markdown** dentro del .txt — solo arte de texto puro.
- Mínimo **40 líneas** salvo que el usuario pida algo mínimo.

## Pasos de ejecución

1. Crear carpeta destino si hace falta.
2. Escribir el `.txt` completo con la herramienta de archivos.
3. En chat: confirmar ruta + invitación a abrir.
4. **No** volcar el mapa en el mensaje del chat.

## Variantes

| Petición del usuario | Forma del mapa |
|---------------------|----------------|
| "mapa del chat" | Árbol con ramas por turno |
| "resumen artístico" | Mandala ASCII central |
| "cronología" | Línea horizontal con hitos |
| "ideas conectadas" | Red de nodos con puentes |

## Si falta contexto

Si el chat es muy corto, amplía con metáforas visuales del tema — pero marca en el pie: `◈ mapa parcial — sesión breve ◈`.