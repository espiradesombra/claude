# ZypyZape / Quijote / Kuramoto
**Víctor Manzanares Alberola · EPSA-UPV Alcoi · 2026**  
`github.com/espiradesombra/claude`

---

## Què és això en 30 segons

Un sistema de **bateria d'inèrcia sense electroquímica** per a aerogeneradors, en tres capes:

| Capa | Nom | Mecanisme | Resposta |
|------|-----|-----------|----------|
| 1 | **ZypyZape** | Transferència cinètica entre rotors via inversor | < 100 ms |
| 2 | **Quijote** | Modulació asimètrica del bobinat trifàsic (balls) + robatori de gravetat (finestra 30-60°) | < 20 ms |
| 3 | **Kuramoto** | Auto-sincronització via bus AC (sense cable addicional) | emergent |

---

## Fitxers clau

```
zypyzape_minigemelo.py   ← simulació digital twin (comencem aquí)
mdc.py                   ← algoritme de factorització MDC (teoria de nombres)
mrauv_goldbach.py        ← model de densitat de primers (Goldbach)
discriminant.py          ← filtre determinista per compostos
sofi_structure.py        ← classificació Sophie Germain
```

---

## Per a enginyers / empreses

**Tres preguntes concretes:**

1. **Compatibilitat:** El sistema no requereix canviar el generador — només el firmware del convertidor. Podeu integrar `zypyzape_minigemelo.py` en el vostre entorn de simulació (MATLAB/Simulink, Python)?

2. **Eficiència:** Heu calculat quant estalviaríeu en filtres actius si useu inèrcia variable `m_i` per absorbir harmònics en lloc de condensadors?

3. **Estabilitat:** Davant d'un ROCOF de 2 Hz/s, quants ms triga el vostre sistema a desconnectar-se? Quijote guanya ~200 ms extra per impuls gravitatori sincronitzat.

---

## Marc teòric (documents 2026)

- `Kuramoto_ZypyZape_VMA_2026.docx` — base matemàtica, swing equation, camp mig via bus AC
- `Inercia_Gravetat_Girs_VMA_2026.docx` — CM electromagnètic vs físic, òrbites en espiral
- `ZypyZape_Quijote_1_3Angular_2026.docx` — finestra 30-60°, 1/3 angular, taula comparativa vs BESS
- `El_Que_Ja_Tenim_VMA_2026.docx` — justificació econòmica i social

---

## Contacte

Víctor Manzanares Alberola  
Investigador independent · EPSA-UPV Alcoi  
`github.com/espiradesombra/claude`

> *"No es que sea fácil. Es que merece la pena."*
