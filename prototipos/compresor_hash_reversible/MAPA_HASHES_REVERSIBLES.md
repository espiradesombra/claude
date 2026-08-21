# Mapa de hashes (¿reversibles?) en el monorepo

| # | Familia | Dónde | ¿Reversible? | Cómo |
|---|---------|-------|--------------|------|
| 1 | **K3 clásico** | `hashtool-work`, `vma-k3/c/k3hash` | **No** (solo fingerprint) | Digest 32-bit; no recupera input |
| 2 | **K3 + Toffoli** | `vma-k3/vma/k3/toffoli_hash.py` | **Parcial** (puerta) | `state ^ (block & prev)` es involutiva en el bit; cadena + banderitas |
| 3 | **Cinemático 1/0** | `hash cinematico…` | **Sí con meta** | Acelera/decelera; guardas `v,t,rachas,pasadas` |
| 4 | **Oceánico / flotación** | `hash cinematico…/29fea3` | **Sí con params** | A/B/C + gravedad/densidades/pasadas |
| 5 | **Espejo degradante** | `compresor y ultima calculadora/01c451` | **Sí con params** | `R=R0·e^{-k I t cosθ}`; paralizas → hash |
| 6 | **Espejo+bloque (compresor HTML)** | `…/708cf0` | **Sí con diccionario/params** | Hash de reflectividad+bloque+clave; reconstrucción busca/dicc |
| 7 | **WiFi / ondas CSI** | mismos HTML + JS seguridad ondas | **Sí si guardas params físicos** | Latencia/Doppler/multicamino → params → hash |
| 8 | **Geométrica / ofuscación** | `geometrica encriptacion+ofuscacion DEMO` | **Cifrado+hash**, no compresor | K3 firma + ofuscación; otro rol |

## Anidación (sí)

Orden recomendado en el prototipo VHC2:

```text
bloque
  → [1] cinemático     (meta trayectoria)
  → [2] espejo         (meta R0,k,I,θ,t)
  → [3] toffoli-k3     (estado + prev)
  → [4] k3-chain       (H_i = K3(H_{i-1} || nest_blob))
```

Cada capa guarda `{family, digest, meta}`.  
La siguiente hashea el **blob anidado** (digest+meta serializado), no solo el texto crudo.

Así el compresor puede:

- deduplicar por digest exterior,
- verificar capa a capa,
- ampliar con oceánico / wifi como capas opcionales.
