---
name: modo-canto
description: >
  Modo canto con notación vocal para interpretación en audio. Formatea letras,
  rimas y respuestas con pausas (,) y énfasis (! / !!). Usar cuando el usuario
  pida cantar, modo canto, ten green bottles, entonación, audio mode, o
  ejecute /modo-canto. El chat TUI no reproduce audio; la notación prepara el
  texto para TTS, lectura en voz alta o Grok con voz.
metadata:
  short-description: "Canto con notación vocal (, ! !!)"
compatibility: Funciona en chat; la interpretación vocal la hace audio/TTS externo.
---

# Modo Canto

Activa el **modo canto**: toda respuesta musical o lírica se entrega en **notación vocal**, no en prosa plana.

## Límite importante (chat vs audio)

- **Grok Build (este TUI)** tiene **dictado por voz** (entrada), no un modo canto de salida integrado.
- El **scrollback de chat** no transmite entonación, ritmo ni silencios con fidelidad.
- Esta skill compensa eso: el texto lleva marcas que un lector humano, TTS avanzado o **Grok con voz** (app/web) puede interpretar.
- En el chat, añade **una línea breve** al inicio: `🎵 modo canto — lee con ritmo; las marcas guían pausas y énfasis`.

## Sistema de notación

| Marca | Significado | Duración / efecto |
|-------|-------------|-------------------|
| `,` | Pausa corta | ~0.3 s, respiración entre frases |
| `,,` | Pausa media | ~0.6 s, tras caídas o giros |
| `,,,` | Silencio dramático | ~1 s, antes del coro o del número |
| `!` | Énfasis suave | Subir tono levemente, acentuar sílaba |
| `!!` | Énfasis fuerte | Pico de energía, casi grito cantado |
| `?` | Entonación interrogativa | Subida al final de la frase |
| `~` | Ligadura / arrastre | Une sílabas sin corte |
| `\|` | Compás / barra | Marca el pulso (opcional en coros) |

### Reglas de escritura

1. **Nunca** entregues la letra sin marcas si el usuario pidió canto.
2. Coloca `!!` en palabras clave: números, acciones (`falls`, `cae`, `quedan`).
3. Usa `,,` o `,,,` después de eventos dramáticos (botella que cae, silencio).
4. Mantén **ritmo regular**: una línea = una frase musical.
5. Si cantas una ronda (Ten Green Bottles), **decrementa el número** y repite el patrón.

## Plantilla: Ten Green Bottles

```
Ten green bottles,, hanging!! on the wall!,
Ten green bottles,, hanging on the wall!,
And if one green bottle accidentally falls,,,
There'll be nine green bottles,, hanging on the wall!
```

Sustituye el número en cada vuelta (nine → eight → … → none).

## Pasos al responder

1. Confirmar modo: `🎵 modo canto activo`.
2. Escribir la letra **solo con notación vocal** (bloque de código o texto monoespaciado).
3. Opcional: línea `Tempo: moderato | Estilo: ronda infantil`.
4. **No** añadir párrafos largos de explicación; máximo 2 líneas de contexto.
5. Si el usuario pide otra canción, aplicar el mismo sistema de marcas.

## Ejemplo mínimo (3 botellas)

```
Three green bottles,, hanging!! on the wall!,
Three green bottles,, hanging on the wall!,
And if one green bottle accidentally falls,,,
There'll be two green bottles,, hanging on the wall!
```

## Referencia

Ver `references/ten-green-bottles.txt` para la letra completa (10 → 0).