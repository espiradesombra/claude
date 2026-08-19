"""
Ciclo: electroiman objeto <-> recorrido (perneo) bajo el mar
=============================================================
Secuencia VMA:
  1) Iman en el OBJETO  -> densidad efectiva sube -> baja solo
  2) En el fondo: rotar objeto, soltar iman
  3) Iman se FIJA / PERNEA al RECORRIDO (otra parte)
  4) Objeto ligero sube (o se reintegra al enjambre)
  5) La "densidad devuelta al recorrido" = la masa del iman ahora pesa en la guia

Contabilidad: m_im * g * Delta_h ganado al bajar con el iman
              = m_im * g * Delta_h que hay que pagar al reponer el iman arriba
              (ya sea subiendo el iman por el recorrido, o rotando la guia).

El electroiman solo es el interruptor (coste ~0 si electropermanente).
"""

from __future__ import annotations

import json
from pathlib import Path

G = 9.81
OUT = Path(__file__).resolve().parent


def cycle_balance(
    m_obj: float = 40.0,
    m_im: float = 8.0,
    delta_h: float = 20.0,  # m de descenso util
    eta_recup_bajada: float = 0.85,  # fraccion de m g h recuperada al bajar (gen/freno)
    eta_subida_obj: float = 0.90,  # si el objeto sube con ayuda mecanica residual
    eta_lift_im: float = 0.88,  # eficiencia al subir el iman por el recorrido
    E_switch_J: float = 2.0,  # pulsos electropermanente (pegar/soltar) por ciclo
    drag_frac: float = 0.05,  # perdidas hidrodinamicas como fraccion de m_tot g h
):
    """
    Un ciclo idealizado (1 iman, 1 objeto, 1 Delta_h).
    Estados de masa:
      A (arriba): iman en objeto
      B (abajo):  iman se pernea al recorrido; objeto suelto
      Reset:      iman vuelve arriba POR EL RECORRIDO (o la guia rota)
    """
    # --- fase 1: bajar con objeto+iman ---
    m_down = m_obj + m_im
    W_grav_down = m_down * G * delta_h  # energia potencial liberada (disponible)
    W_drag_down = drag_frac * W_grav_down
    W_rec_down = eta_recup_bajada * (W_grav_down - W_drag_down)  # lo que sacas a bus/almacen

    # --- fase 2: soltar iman abajo (switch) ---
    E_switch = E_switch_J

    # --- fase 3: objeto solo sube (flotabilidad o enjambre lo reintegra) ---
    # Si flota "gratis" por Arquimedes, NO recuperamos m_obj*g*h al subir:
    # la flotabilidad hace el trabajo (el fluido "paga" y al comprimir/empujar se reequilibra
    # en el gran sistema oceano; para el modulo local el subida no nos cobra motor).
    # Pero tampoco podemos CONTAR m_obj*g*h como ganancia neta extra: ya la descontamos
    # en el potencial del sistema completo.
    # Contabilidad del MODULO (caja):
    #   - Al bajar con m_down recuperamos parte de m_down*g*h
    #   - Al subir el objeto ligero sin motor: 0 coste electrico en el objeto
    #   - El iman se quedo abajo: el RECORRIDO debe subir m_im
    W_lift_im_ideal = m_im * G * delta_h
    W_lift_im_elec = W_lift_im_ideal / eta_lift_im  # motor de guia / rotacion recorrido

    # Opcional: si el objeto NO flota solo y hay que remolcarlo arriba
    W_tow_obj = 0.0  # 0 = flota solo / enjambre lo recoge sin coste contabilizado aqui

    # Energia electrica neta del ciclo (positivo = generamos, negativo = consumimos)
    E_out = W_rec_down
    E_in = E_switch + W_lift_im_elec + W_tow_obj
    surplus = E_out - E_in

    # Desglose "ilusion":
    # Si alguien olvida subir el iman: parece que gana m_im*g*h
    surplus_if_forget_im_reset = E_out - E_switch  # TRAMPA contable

    return {
        "m_obj": m_obj,
        "m_im": m_im,
        "delta_h_m": delta_h,
        "W_grav_down_J": W_grav_down,
        "W_rec_down_J": W_rec_down,
        "W_lift_im_ideal_J": W_lift_im_ideal,
        "W_lift_im_elec_J": W_lift_im_elec,
        "E_switch_J": E_switch,
        "E_out_J": E_out,
        "E_in_J": E_in,
        "surplus_J": surplus,
        "surplus_if_forget_im_reset_J": surplus_if_forget_im_reset,
        "is_generator": surplus > 0,
        "note": (
            "Pernear el iman al recorrido SOLO mueve el peso al carril. "
            "Quien cierre el ciclo (rotar guia / subir iman por la otra parte) "
            "paga m_im*g*Dh. Electroiman ~ interruptor barato."
        ),
    }


def main():
    print("=" * 66)
    print("IMAN EN OBJETO -> BAJA -> PERNEA AL RECORRIDO (otra parte)")
    print("=" * 66)

    base = cycle_balance()
    print("\nCaso base (m_obj=40, m_im=8, Dh=20 m, eta_bajada=0.85, eta_lift_im=0.88):")
    for k in [
        "W_grav_down_J",
        "W_rec_down_J",
        "W_lift_im_ideal_J",
        "W_lift_im_elec_J",
        "E_switch_J",
        "surplus_J",
        "surplus_if_forget_im_reset_J",
        "is_generator",
    ]:
        v = base[k]
        print(f"  {k}: {v:.2f}" if isinstance(v, float) else f"  {k}: {v}")

    print("\nBarrido eta_lift_im (eficiencia de devolver el iman arriba por el recorrido):")
    print(f"  {'eta_lift':>8} {'surplus_J':>12} {'gen?':>6}")
    rows = []
    for eta in [0.50, 0.70, 0.85, 0.92, 0.98, 1.00]:
        r = cycle_balance(eta_lift_im=eta)
        rows.append(r)
        print(f"  {eta:8.2f} {r['surplus_J']:12.1f} {'SI' if r['is_generator'] else 'NO':>6}")

    print("\nBarrido m_im (mas iman = mas 'ayuda' a bajar = mas peaje de subida):")
    print(f"  {'m_im':>6} {'W_rec':>10} {'W_lift_im':>10} {'surplus':>10}")
    for mi in [2.0, 5.0, 8.0, 15.0, 25.0]:
        r = cycle_balance(m_im=mi)
        print(
            f"  {mi:6.1f} {r['W_rec_down_J']:10.1f} {r['W_lift_im_elec_J']:10.1f} "
            f"{r['surplus_J']:10.1f}"
        )

    # Caso ideal sin perdidas: surplus debe ser ~0 o negativo por drag
    ideal = cycle_balance(
        eta_recup_bajada=1.0, eta_lift_im=1.0, E_switch_J=0.0, drag_frac=0.0
    )
    print("\nIdeal (eta=1, sin drag, switch=0):")
    print(f"  surplus = {ideal['surplus_J']:.4f} J  (debe ser ~0)")
    print(f"  W_rec_down = {ideal['W_rec_down_J']:.2f}")
    print(f"  W_lift_im  = {ideal['W_lift_im_elec_J']:.2f}")
    print(
        f"  Nota: W_rec incluye m_obj*g*h + m_im*g*h; "
        f"solo se devuelve m_im*g*h al subir el iman."
    )
    print(
        "  El m_obj*g*h 'recuperado' al bajar se 'gasta' conceptualmente cuando el\n"
        "  objeto flota/sube: la boyancia hace trabajo igual a m_obj*g*h en el fluido.\n"
        "  Si contamos SOLO la caja electrica del modulo y asumimos subida de objeto\n"
        "  gratis por flotacion, HAY que restar m_obj*g*h del saldo (o no contarlo\n"
        "  como E_out). Correccion contable abajo."
    )

    # Contabilidad CORRECTA de caja electrica + flotacion
    def cycle_balance_honest(**kw):
        r = cycle_balance(**kw)
        # Al bajar recuperamos (m_obj+m_im)gh * eta
        # Al flotar el objeto, no pagamos motor, pero NO es energia electrica gratis
        # reutilizable sin limite: el potencial del objeto se restauro via boyancia.
        # Para no doble-contar: E_out_util solo la parte atribuible al iman que
        # luego pagamos al subir, mas cualquier exceso real de fuentes externas.
        # Honest split:
        m_obj, m_im, dh = r["m_obj"], r["m_im"], r["delta_h_m"]
        eta_b = kw.get("eta_recup_bajada", 0.85)
        drag = kw.get("drag_frac", 0.05)
        W_obj = m_obj * G * dh
        W_im = m_im * G * dh
        W_drag = drag * (W_obj + W_im)
        # recuperacion repartida proporcional a masas
        W_rec_obj = eta_b * (W_obj - drag * W_obj)
        W_rec_im = eta_b * (W_im - drag * W_im)
        # la parte obj se anula con la subida flotante (no es surplus electrico durable)
        # la parte im se anula con lift del iman (con eta_lift)
        eta_l = kw.get("eta_lift_im", 0.88)
        E_switch = kw.get("E_switch_J", 2.0)
        surplus_honest = W_rec_im - (W_im / eta_l) - E_switch
        r["W_rec_obj_J"] = W_rec_obj
        r["W_rec_im_J"] = W_rec_im
        r["surplus_honest_J"] = surplus_honest
        r["is_generator_honest"] = surplus_honest > 0
        return r

    print("\n--- Contabilidad HONESTA (no contar m_obj*g*h como generacion neta) ---")
    print(f"  {'eta_lift':>8} {'rec_im':>10} {'lift_im':>10} {'surplus_h':>10}")
    honest_rows = []
    for eta in [0.70, 0.85, 0.92, 0.98, 1.00]:
        r = cycle_balance_honest(eta_lift_im=eta)
        honest_rows.append(r)
        print(
            f"  {eta:8.2f} {r['W_rec_im_J']:10.1f} {r['W_lift_im_elec_J']:10.1f} "
            f"{r['surplus_honest_J']:10.1f}"
        )

    h1 = cycle_balance_honest()
    print(f"\n  Caso base surplus_honest = {h1['surplus_honest_J']:.1f} J -> generador? {h1['is_generator_honest']}")

    print(
        """
LECTURA DE TU SECUENCIA
-----------------------
  [arriba]  electroiman EN EL OBJETO  ->  mas denso  ->  baja solo
  [abajo]   rotas objeto, sueltas iman  ->  iman PERNEA AL RECORRIDO (otra parte)
  [objeto]  queda casi neutro / flotante  ->  sube o se reintegra al enjambre
  [recorrido] ahora CARGA con la masa del iman en la cota baja

  "Devolver la densidad al recorrido" = el peso del iman pasa del objeto a la guia.
  No desaparece. La otra parte del recorrido (o su rotacion) es quien debe
  reponer el iman arriba. Ahi se paga m_im * g * Delta_h.

  Por eso:
    - Soltar/pegar con electropermanente ~ barato (interruptor)
    - Perneo al recorrido = embrague mecanico util (desperneado / fase)
    - NO es robo de gravedad
    - SI es interesante como logistica de masa y buffer de densidad
"""
    )

    out = {
        "base_naive": base,
        "base_honest": h1,
        "verdict": {
            "roba_gravedad": False,
            "electroiman_es_fuente": False,
            "electroiman_es_interruptor": True,
            "perneo_al_recorrido_mueve_peso_a_la_guia": True,
            "quien_paga": "subir/reponer el iman en la cota alta (recorrido, rotacion o enjambre)",
            "utilidad": "acoplamiento, densidad conmutable, perneo/desperneo, no generador primario",
        },
    }
    path = OUT / "ciclo_iman_perneo_recorrido.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"JSON: {path}")


if __name__ == "__main__":
    main()
