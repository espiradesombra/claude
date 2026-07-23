# Informe de verificació + inicios Quijote / ZypyZape
Data: 2026-07-23
Pack: VMA_mates_rescat_2026

## 1. Verificació de codi (Python 3.12)

| Mòdul | Estat | Detall |
|-------|-------|--------|
| cribas.py Desmemoriada | OK | prims ≤200 exactes (46) |
| cribas.py Híbrida | OK | prims ≤200 exactes |
| cribas.py Modular6k | OK (fix 2026-07-23) | fase 2p/4p corregida; match exacte ≤500 |
| criba_hibrida.py mini | OK | prims ≤100 exactes |
| criba_modular.py mini | PARCIAL | marca candidats; no és llista neta de prims |
| criva.py | OK | criva(1000) dens≈0.089 |
| mdc.py mdc_factor | OK | 91→(7,13); 143→(11,13); 1147→(31,37) MDC hit |
| mdc_v23.py | OK | carrega (check_factor, bateries…) |
| newton_rapido.py | PARCIAL | demo ok; convergència j per (E=100,b=10) dèbil en test ràpid |
| fermat_modular.py | OK | identitats F0…Fn = F(n+1)-2 ✓; alineació modular OK |
| sofi_structure.py | OK | classify_sofi(30) OK |
| estructura_sofi mini | OK | LSG [5,11,23,29] |
| siguiente_primo.py | OK (fix 2026-07-23) | roda 6k±1 per defecte; karnaugh experimental opcional |
| gemelo_v942 / v10 | OK syntax | compile OK (cal matplotlib per executar gràfics) |
| anexoF OpenMP .c | NO PROVAT | gcc no està al PATH d’aquest PC |

### Resum (post-fix)
- **Tot el nucli Python de cribes + siguiente_primo passa validació.**
- **XFI gemelo:** `10_inicio_quijote_zypyzape/gemelo_xfi_avion_4turbinas.py` (+ PNG).
- **C:** instal·lar MinGW/gcc o WSL si vols OpenMP.

### XFI — avió 4 turbinas
```
python 10_inicio_quijote_zypyzape\gemelo_xfi_avion_4turbinas.py
→ xfi_avion_dinamica.png
```
Model conceptual (ondulació climb/dive, rols gen/thr/buf). No és perpetuum mobile.

## 2. Inici QUIJOTE = «molinos en molinos»

No hi ha un fitxer amb el títol exacte «molinos en molinos», però el concepte és clar i documentat:

**Quijote = molí dins del molí:** masses desplaçables a les pales del rotor eòlic
que modulen la inèrcia J(r) = J_G + N_b · m_q · r².

### Fitxers d’inici (canònics)
- `10_inicio_quijote_zypyzape/01_CONTEXT_ZYPYZAPE.txt` — secció QUIJOTE
- `10_inicio_quijote_zypyzape/02_MATH_QUIJOTE_3vs7.txt` — 3 pales vs 7, J constant
- Codi madur: `09_gemelos_zypyzape_ultims/zypyzape_twin_v4_8_quijote.py`
- Codi paper: `09_.../gemelo_v94.py` … `gemelo_v942.py`, `gemelo_v10_cp_dinamic.py`
- Carpeta original completa: `monton viejo revisar\archivos\claude-main-extract\claude-main\Quijote\`

Física clau:
- Cap a fora → carrega bateria (J↑)
- Cap a dins → descarrega
- Equació: J·dω/dt = T_net − ω·(dJ/dt)

## 3. Inici ZYPYZAPE = «avió eterno / 4 turbinas»

### A) Idea conceptual «avió / 4 molinos» (xat, no codi formal)
Origen: `txt25 5\Ver nuevas publicaciones.txt`  
Extract al pack: `10_inicio_quijote_zypyzape/INICIO_avion_4turbinas_extract.txt`

Idees originals de Víctor (parafrasejades del xat):
1. **Cotxe/camió amb 4 molinos** + ZypyZape com a bateria d’inèrcia  
   (v_coche + v_acumulada + v_viento; Lenz vs termodinàmica).
2. **Avió amb 3 motors** connectats a l’estil ZypyZape:  
   en baixada capta cinètica; en picat accelera; ales tornen energia en alçada  
   → «avió que no reposta» (drones/càmeres, experimental).

### B) Codificació formal ZypyZape (evolució)
No hi ha un `avion_eterno_4turbinas.py` dedicat. L’enginyeria formal va anar a **parcs eòlics**:

| Versió | Topologia | Fitxer |
|--------|-----------|--------|
| INICI twin | N=5 (parell central + anell 3) | `zypy_zape_digital_twin_INICIO.py` |
| v2–v3 | 5 nodes document tècnic | `zypyzape_twin_v2.py` … |
| minigemelo | N_TURB=5 | `09_.../zypyzape_minigemelo.py` |
| v4.8 + Quijote | mòdul 10 turbines (2+2×4), N_q=5 amb Quijote | `zypyzape_twin_v4_8_quijote.py` |
| Papers gemelo | Quijote+ZZ v9.4.x / v10 | `gemelo_v94*.py` |

Paràmetres tipus (context v4.8):
- 2.5 MW/turbina, R=60 m, J=5e6 kg·m², cicle 0.4 Hz
- Topologia: parelles que s’intercanvien energia (una accelera, l’altra frena)

### C) «1–4 turbinas moduladores»
En docs de parc 10×5 MW: grup òptim 6–9 en MPPT + **grup modulador 1–4** precarregades  
(extracte Libro4 / Apiñon). Això és el lligam més proper a «4 turbinas» en codi/docs tècnics.

## 4. On començar a llegir (ordre recomanat)

1. `10_inicio_quijote_zypyzape/01_CONTEXT_ZYPYZAPE.txt`
2. `10_inicio_quijote_zypyzape/02_MATH_QUIJOTE_3vs7.txt`
3. `10_inicio_quijote_zypyzape/INICIO_avion_4turbinas_extract.txt` (idea avió/cotxe)
4. `10_.../zypy_zape_digital_twin_INICIO.py` (primer twin N=5)
5. `09_.../zypyzape_twin_v4_8_quijote.py` (estat més madur)
6. `09_.../gemelo_v942.py` (paper consolidat)

## 5. Properes accions opcionals
- Fix CribaModular6k + siguiente_primo
- Portar idea «avió 3/4 motors» a un gemelo nou a part del parc eòlic
- Instal·lar gcc per provar anexoF OpenMP
