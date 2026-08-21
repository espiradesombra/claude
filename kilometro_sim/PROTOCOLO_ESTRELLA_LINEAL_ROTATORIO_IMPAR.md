# Protocolo estrella — lineal + rotatorio + rearme impar

**Cómo lo describe VMA (sin aplastar a “solo batería”):**

> Es el final del *hurto gravitatorio*: **se puede**, pero es **complicado de narices**.  
> Media bajada **lineal** → media **rotatoria** del recorrido → ir **robando distancia** a ese ciclo → **rearme cada N impar** de vueltas del recorrido por **flotación + perneo de pesos**.  
> Rol: **estrella de renovables** = buffer + convertidor de fase que saca **generador además del inventario**, si el faseo impar cierra bien.

No confundir con `sim_enjambre_pesos_impar` (ahí “impar” = tamaño de stock).  
Aquí **impar** = ritmo de **rearme respecto al recorrido que gira**.

---

## 1. Geometría del ciclo (una vuelta lógica)

```text
        ESTACIÓN S (siempre el mismo lado)
        recarga distancia / perneo
                 │
    ┌────────────┴────────────┐
    │                         │
    ▼                         ▼
 MEDIA LINEAL              MEDIA ROTATORIA
  bajada Δh                 recorrido gira φ
  lastre enganchado         sinfín / eje
  → GENERADOR               → roba / rearma DISTANCIA
  (PE lastre → julios)        (geometría del Kilómetro)
    │                         │
    └────────────┬────────────┘
                 │
        cada N = 1,3,5… vueltas:
        FLOTACIÓN + PERNEO
        (objeto ligero sube; pesos cambian de extremo)
```

| Tramo | Qué pasa | Canal energético |
|-------|----------|------------------|
| **Bajada media lineal** | Objeto + lastre extra descienden en guía casi recta | **Generador** (freno regenerativo / sinfín) |
| **Media rotatoria** | El recorrido gira; se “roba distancia” (rearmar radio/fase sin pagar toda la PE otra vez) | **Peaje de actuador / inercia** — no es el almacén |
| **Rearme impar** | Cada N impar: flotar (n≤3), pernear pesos al otro extremo / stock | **Batería de inventario** + reset del objeto |

---

## 2. Dos cuentas (obligatorio no mezclarlas)

| Cuenta | Qué mide | Cuándo crece |
|--------|----------|--------------|
| **Batería** | Lastres que bajan de cota ALTA→BAJA (o no se restituyen) | Cada patada que aparca kg abajo |
| **Generador** | Trabajo eléctrico del eje/sinfín en la media lineal (y lo regenerable en rotatorio) | Mientras hay bajada / fase favorable |
| **De más** | Generador que **no** se explica solo por vaciar batería, a **SOC casi fijo** tras rearmes impares | La tesis estrella |

Si solo miras julios mientras vacías el desván → estás midiendo batería.  
Si tras rearmes impares el inventario vuelve y aún hay neto eléctrico útil → eso sería el **plus**.

---

## 3. Por qué “complicado de narices”

1. Sincronizar **fase lineal** con **φ del tubo** (no adelantar ni retrasar el perneo).  
2. Make-before-break de pernos bajo carga hidrodinámica.  
3. Calibrar **3 flota / 4 se hunde** con tolerancia real de densidad.  
4. Robar distancia **sin** gastar en actuador más de lo regenerado en lineal.  
5. Rearme **solo en impares**: el patrón par/impar respecto a cos(φ) no es simétrico — hay que **medir**, no intuir.  
6. Cerrar 1.ª ley en laboratorio: dos vatímetros + celdas de carga en lastres.

Viable como ingeniería dura (TRL bajo→medio).  
No es un PDF de sobreunidad: es un **ciclo híbrido** que hay que instrumentar.

---

## 4. Scripts

| Archivo | Rol |
|---------|-----|
| `sim_recorrido_rota_patadas_impar.py` | Barrido N par/impar; modo A (almacén) vs B (SOC fijo) |
| `sim_enjambre_pesos_impar.py` | Solo inventario 3/4 (no es este protocolo) |
| `doble_km_90_perneo.py` | Buffer de picos / perneo inercial |
| Este markdown | Especificación del protocolo estrella |

### Cómo leer el experimento actual

```bash
python kilometro_sim/sim_recorrido_rota_patadas_impar.py
```

- **Modo A**: puede lucir julios por vaciar batería.  
- **Modo B (SOC fijo)**: prueba dura del “generador de más”.  
  En el modelo 1.ª-ley + lift de la v1 del script, el plus **no aparece** (gen_de_mas &lt; 0).  
  Eso **no cierra** tu diseño: el script aún **no** separa bien “media lineal vs media rotatoria” ni el rearme impar por flotación como máquina de estados.  
  Es el **esqueleto de contabilidad**; el siguiente paso de código es la máquina de estados de §1.

---

## 5. Variante: rearme cada 3 ciclos + robar 1/4 de recorrido

**Opinión de diseño (VMA):**

| Idea | Lectura |
|------|---------|
| **Pernar cada 3 ciclos** | Cadencia más lenta que “cada impar”: deja acumular 2 ciclos de generación “abiertos” y un 3.º de reset. Buen ritmo de enjambre (menos conmutaciones → menos `E_perno`). Riesgo: stock BAJA crece 2 ticks antes del reset. |
| **Robar 1/4 de recorrido por patada** para **potencial de flotación** | En la media rotatoria, cada patada adelanta/retiene ~90° de fase o ~1/4 de longitud efectiva, de modo que al soltar lastre el módulo quede en zona donde la boyancia hace más trabajo (subida barata). Es **faseo hidrostático**, no julios gratis: estás eligiendo *cuándo* el neutro actúa. |

Encaja con el gemelo nuevo `ARRIBA Red Renovable/…reset-de-potencial-de-flotación…`: el reset no es lift bruto siempre; es **perneo + geometría** para recuperar potencial de flotación.

**Criterio de éxito:** bajar `E_rob_dist` por patada (objetivo ~1/4 de peaje de una media rotatoria completa) sin perder `W_gen` de la lineal.

---

## 6. Frase de posicionamiento (renovables)

> Kilómetro estrella = **mitad lineal que genera**, **mitad rotatoria que rearma distancia**, **rearme cada 3 (o impar) por flotación y perneo**.  
> Batería de lastre + convertidor de fase. El hurto gravitatorio “se puede” como **ciclo híbrido bien faseado**, no como pozo infinito.

---

*Documento de diseño VMA / kilometro_sim, 2026-08.*
