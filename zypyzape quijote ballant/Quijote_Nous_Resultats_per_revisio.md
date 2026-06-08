# Quijote — Resultats nous per a revisió crítica

**Autor:** Víctor Manzanares Alberola — EPSA, UPV (Alcoi)
**Repositori:** github.com/espiradesombra/claude
**Estat del document:** esborrany per a revisió externa (GPT i altres)

---

## Petició explícita al revisor

Vull una **revisió dura**, no una validació. Concretament:

1. Busca errors de física, de signe, de doble comptabilització o de
   formulació abans que resultats positius.
2. Qüestiona els supòsits del model (els marco com a supòsits, no com a
   mesures).
3. Digues-me quines conclusions són **estructurals** (independents del
   soroll numèric) i quines depenen de la calibració concreta.
4. Si alguna cosa es pot demostrar analíticament, prefereixo la demostració
   a la simulació.

Hi ha dos blocs nous: (A) la comparativa d'estabilitat de freqüència contra
l'estat de l'art, i (B) l'estudi del bombeig paramètric aplicat al ball de
pesos amb el pes a la punta. Tots dos parteixen d'un fet ja tancat en treball
anterior: **l'hurto gravitatori en cicle tancat és nul per conservació de
l'energia** (demostrat per tres camins: funció d'estat, Lagrange, ∮F·dr=0,
amb residu de balanç que convergeix a zero en refinar la malla).

---

## BLOC A — Comparativa d'estabilitat de freqüència

### Objectiu

Situar Quijote i ZypyZape **no en energia** (ja sabem que no hi ha guany
energètic estructural) sinó en **dinàmica de xarxa**, contra l'estat de l'art
real: control electrònic amb inèrcia sintètica i bateries grid-forming.

### Tecnologies comparades

Parc de 44 turbines NREL 5MW, mateix esdeveniment de pèrdua de generació:

1. **BASE** — MPPT pur, sense suport de freqüència.
2. **Inèrcia sintètica** — suport electrònic via convertidor, resposta en
   mil·lisegons (estat de l'art).
3. **Quijote** — massa mòbil a la punta, resposta mecànica de ~0,3–3 s.
4. **ZypyZape** — coordinació de fase tipus Kuramoto entre turbines.
5. **ZypyZape + Quijote** — combinació.
6. **Bateria grid-forming (BESS)** — referència de mercat, resposta en ms.

### Mètriques (jutge)

RoCoF (taxa de canvi de freqüència, 0–500 ms), nadir de freqüència, energia
aportada durant l'esdeveniment, i dispersió de velocitats angulars del parc
(sincronisme).

### Resultat honest i limitació del model

**Important — limitació admesa:** el model de freqüència de xarxa que vaig
usar (swing equation agregada amb amortiment) no està prou ben calibrat per
discriminar numèricament les tecnologies: amb la inèrcia base assignada, totes
donen un RoCoF pràcticament idèntic i nadirs que difereixen només a la quarta
xifra. **No considere fiables els valors numèrics absoluts d'aquesta taula.**
Per a una taula numèrica fiable caldria una eina d'enginyeria de xarxa
professional (PSS/E, DIgSILENT, o un swing multi-màquina calibrat amb dades
reals d'un TSO).

El que **sí** és estructural, i no depèn d'aquesta calibració, és la
separació d'escales de temps, que es pot afirmar des de la física:

| Mètrica | Qui guanya | Raó (estructural) |
|---|---|---|
| RoCoF (0–500 ms) | Electrònica (BESS, synth-inertia) | Responen en ms; Quijote mecànic respon en segons. Separació d'escales temporals. |
| Nadir (banda de segons) | Electrònica, amb marge | La caiguda ja queda "pre-aplanada" abans que Quijote desplegue la massa. |
| Energia aportada | Cap (net) | Conservació: el suport és transferència, no generació. |
| Sincronisme de parc (dispersió ω) | ZypyZape | És coordinació (Kuramoto), no velocitat. Única mètrica amb valor diferencial. |

### Conclusió del bloc A

En estabilitat de freqüència, l'electrònica de potència guanya en les
mètriques que es paguen (RoCoF, nadir) perquè actua en una escala temporal
inferior a la del problema. El valor diferencial de **ZypyZape** és la
**coordinació de parc** (sincronisme), perquè és software i no depèn de la
velocitat mecànica. **Quijote** no aporta avantatge competitiu en cap mètrica
clau.

Aquest resultat encaixa amb el límit que proposes formalitzar:

> τ_actuador ≫ τ_grid  ⟹  contribució del mecanisme mecànic → 0

És a dir: quan el sistema base ja té control electrònic ràpid, l'espai
funcional de la massa mecànica desapareix. No perquè Quijote "no funcione",
sinó perquè el problema ja està resolt per un mecanisme més ràpid i barat.

**Pregunta oberta al revisor:** et sembla correcte afirmar la taula
"qui guanya" només des de l'argument d'escales temporals, sense dependre de la
simulació mal calibrada? O creus que cal el model professional per sostindre-la?

---

## BLOC B — Bombeig paramètric al ball de pesos (pes a la punta)

### Motivació: el columpi

Un xiquet es gronxa sol mitjançant **bombeig paramètric**: puja i baixa el seu
centre de massa dues vegades per oscil·lació, sincronitzat amb la fase. Açò
amplifica l'oscil·lació de veritat. La pregunta natural: **es pot aplicar el
mateix al ball de pesos de Quijote per extraure energia?**

La hipòtesi a comprovar: amb el pes a la punta, una massa del 8% de la pala, i
moviments molt xicotets sincronitzats (bombeig a 2× la freqüència orbital),
l'hurto és **positiu, negatiu o nul**?

### Paràmetres

- Pes a la punta: R = 54 m.
- Massa per pala: 8% de la pala NREL (≈ 1.440 kg, sobre pala de ~18.000 kg).
- Moviment xicotet: Δr = 0,05 m (5 cm).
- Modulació paramètrica: r(t) = R − Δr·sin(2θ + φ), escombrant la fase φ.
- ω **lliure** (integració ODE completa de Lagrange amb DOP853).
- Balanç energètic complet: ΔE_rotació, treball d'actuador, drift de ω.

### Resultat

Escombrant 8 fases de bombeig (buscant el cas més favorable):

- En **totes** les fases, ΔE_rotació i W_actuador són **iguals al cèntim**.
  Quan la rotació guanya +0,56 kJ, és perquè l'actuador hi ha posat +0,56 kJ.
  Quan perd, és perquè l'actuador ho ha recuperat.
- El **drift de ω és 0,00%** en tots els casos: el cicle tanca de veritat
  (la velocitat angular torna al valor inicial).
- **Hurto net = 0** en totes les fases provades.

### Interpretació (la clau)

El bombeig paramètric **funciona** com a mecanisme de transferència (per això
ω es mou de veritat segons la fase), exactament com el columpi. Però l'energia
ix de l'**actuador**, no de la gravetat — igual que en el columpi ix del
**múscul** del xiquet, no de la gravetat. La gravetat actua com a mediadora,
mai com a font.

La distinció fina amb el columpi: el gronxador acumula energia perquè és una
oscil·lació **lliure que creix** (el cicle d'energia no es tanca; el xiquet es
gronxa cada vegada més alt sense tornar l'energia). Quijote en règim
estacionari **tanca el cicle cada volta** (la massa torna, ω torna), i per
funció d'estat ΔE_gravitatori = 0 exacte.

### Conclusió del bloc B

L'hurto és **NUL**, no positiu ni negatiu. El bombeig paramètric d'inèrcia
variable en un rotor no extrau energia gravitatòria neta en règim estacionari;
només media una transferència actuador ↔ rotació. Aquesta és l'última via
sofisticada per intentar l'hurto (la que fan els columpis de veritat), i també
dóna zero, per la mateixa raó fonamental que totes les anteriors
(commutació radial, camarón, velocitat asimètrica, pes a la punta).

**Possible resultat publicable (honest):**

> "El bombeig paramètric d'inèrcia variable en un rotor no extrau energia
> gravitatòria neta en règim estacionari (cicle tancat); només media
> transferència d'energia entre l'actuador i la rotació, anàlogament al
> bombeig muscular d'un columpi. Verificat amb integració dinàmica de ω
> lliure i balanç energètic tancat (drift de ω < 0,01%)."

**Preguntes obertes al revisor:**

1. La distinció "oscil·lació lliure que creix (columpi)" vs "règim
   estacionari de cicle tancat (Quijote)" — la trobes ben plantejada, o hi ha
   un règim intermedi (rotor parcialment lliure) que caldria explorar?
2. El càlcul del treball d'actuador usa una aproximació de "pala mitjana"
   (Q_r per pala promitjat). Creus que açò amaga algun terme creuat entre
   pales que caldria tractar pala a pala?
3. Hi ha alguna configuració de bombeig (fase, freqüència múltiple, no
   sinusoïdal) que escape de la nul·litat, o la conservació la prohibeix
   sense excepció (com crec)?

---

## Resum global per al revisor

Després d'explorar **totes** les vies conegudes — commutació radial segons
fase, camarón (encollir per accelerar), velocitat asimètrica, pes a la punta
amb moviments curts, i bombeig paramètric — el resultat és consistent:

- **Hurto gravitatori:** nul, per conservació (demostrat per múltiples camins).
- **Captació de vent extra:** efecte mitjà ≈ 0, dominat pel soroll (μ ≪ σ).
- **Suport de freqüència vs estat de l'art:** dominat per l'electrònica;
  Quijote arriba tard (separació d'escales temporals).
- **Valor diferencial real:** només ZypyZape com a **coordinació de parc**
  (sincronisme), que és software i no depèn de velocitat mecànica.

El treball, per tant, canvia de categoria: no és tecnologia de generació ni
millora d'eficiència, sinó un **estudi rigorós de control multi-escala i del
límit físic de l'actuació retardada** en estabilitat de xarxa, més una
**refutació neta i ben documentada** d'una hipòtesi de tipus moviment perpetu.

Agrairé que ataques especialment els punts marcats com a limitacions
(calibració del model de xarxa al bloc A; aproximació de pala mitjana al
bloc B) i que em digues si les conclusions estructurals se sostenen sense
dependre d'aquests punts febles.
