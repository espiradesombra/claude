# VMA — reglas del proyecto

- Autor: Víctor Manzanares Alberola (VMA).
- Motor principal: `K3CrossVerifier` en `vma/k3/cross_verifier.py`.
- Firma anti-tamper: `K3_FIRMA_LOGICA` siempre como bloque cero.
- Tres modos: Usual (32/0), Propio A (16/stride), Propio B (marcas finales).
- K3 NO es criptográfico — integridad industrial / telemetría, no contraseñas.
- Código C legado en `c/k3hash/` — librería básica streaming, sin firma ni stride.
- Preferir Python para prototipos; C para embebido y binarios finales.
- Nombres de proyecto: Theorem K3, Proyecto 33x1, Gap Prime, Fase, Banderita.