"""
Simulacion: enjambre de pesos + objeto 3-flota / 4-se-hunde
==========================================================
- Recarga de distancia SIEMPRE en el mismo lado (estacion S, cota alta).
- Cada peso: 2 pernos (R recorrido / O objeto). Make-before-break.
- Stock finito de pesos en guia ALTA. Al bajar, el peso se aparca ABAJO.
- N = 1,2,3,4,... ciclos: se ve el "impar XD" = inventario, no sobreunidad.

Salidas:
  - tabla por ciclo y acumulados
  - PNG de energia vs stock
  - JSON resultados
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent
G = 9.81


@dataclass
class Params:
    m_obj_base: float = 30.0  # kg  (cuenta como "masa del cuerpo"; flotacion calibrada aparte)
    m_peso: float = 10.0  # kg cada peso discreto
    # Con n_pesos en el objeto:
    #   n <= n_float  -> flota
    #   n >= n_sink   -> se hunde
    n_float: int = 3
    n_sink: int = 4
    # Pesos fijos en las puntas del objeto (siempre van con el objeto)
    n_puntas: int = 2
    # Stock inicial en guia ALTA (disponibles para enganchar el 4.o, 5.o...)
    n_stock_alta: int = 5
    n_stock_baja: int = 0
    delta_h: float = 15.0  # m  cota alta -> baja
    eta_gen: float = 0.85  # recuperacion electrica en bajada
    eta_lift: float = 0.90  # si reponemos pesos subiendo (modo opcional)
    E_perno_J: float = 1.5  # coste electrico por conmutacion de un peso (2 pernos)
    drag_frac: float = 0.06  # perdidas hidrodinamicas fraccion de m g h
    # n_pesos en objeto al inicio = n_puntas + extras enganchados (0)
    # Para flotar al inicio necesitamos n_puntas <= n_float y sin extra
    # "con 3 flota con 4 se hunde": interpretamos n_total_pesos_en_objeto


@dataclass
class State:
    n_alta: int
    n_baja: int
    n_on_obj: int  # pesos en el objeto (incluye puntas)
    h_obj: str = "alta"  # "alta" | "baja"
    E_out: float = 0.0  # energia electrica recuperada acumulada
    E_in: float = 0.0  # energia electrica gastada (pernos + lifts opcionales)
    cycles_done: int = 0
    log: list = field(default_factory=list)

    @property
    def surplus(self) -> float:
        return self.E_out - self.E_in

    def floats(self, p: Params) -> bool:
        return self.n_on_obj <= p.n_float

    def sinks(self, p: Params) -> bool:
        return self.n_on_obj >= p.n_sink


def mgh(m: float, dh: float) -> float:
    return m * G * dh


def run_cycle(state: State, p: Params, do_lift_reset: bool = False) -> bool:
    """
    Un ciclo de 'recarga de distancia' en lado S (siempre alta):
      1) En cota alta, enganchar 1 peso de stock_alta al objeto (si hace falta para hundir)
      2) Bajar generando
      3) En cota baja, pernear el peso extra al recorrido (stock_baja++)
      4) Objeto con n_puntas (o n_float) sube 'gratis' por flotacion a cota alta
      5) Recarga de distancia en el MISMO lado S (geometria) — coste ~0 contable aqui
      6) Opcional: subir un peso de baja->alta (cierra almacen)

    Devuelve False si no puede completar el ciclo (sin stock alta y no hunde).
    """
    events = []

    # --- 0) objeto debe estar arriba y flotar ---
    if state.h_obj != "alta":
        events.append("FAIL: objeto no esta en cota alta")
        state.log.append({"cycle": state.cycles_done + 1, "ok": False, "events": events})
        return False

    # asegurar solo pesos de puntas al inicio de ciclo (estado "3 flota")
    # extras se enganchan ahora
    extras = state.n_on_obj - p.n_puntas
    if extras < 0:
        state.n_on_obj = p.n_puntas
        extras = 0

    # --- 1) enganchar pesos hasta poder hundirse (tipicamente 1 extra: 2+1+? -> 4) ---
    # Necesitamos n_on_obj >= n_sink. Partimos de n_puntas (2) + ya enganchados.
    # Usuario: con 3 flota, con 4 se hunde.
    # Modelo: n_on_obj cuenta TODOS los pesos en el objeto.
    # Inicio de ciclo util: n_on_obj = 3 (flota) o n_puntas=2 y un "peso de cuerpo" virtual.
    # Para cuadrar 3/4: n_on_obj inicial = 3 (dos puntas + 1 fijo de chasis contado como peso).
    # Mas simple: n_min_float = 3 fijos en objeto; el 4.o se coge del stock.
    # Reinterpretamos: n_on_obj al reposo = p.n_float (3), y cogemos 1 del stock para hundir.
    if state.n_on_obj < p.n_float:
        state.n_on_obj = p.n_float  # calibracion de reposo flotante

    need = p.n_sink - state.n_on_obj  # suele ser 1
    if need < 0:
        need = 0

    if state.n_alta < need:
        events.append(
            f"STOP: stock alta={state.n_alta} < need={need} para hundir "
            f"(n_on_obj={state.n_on_obj})"
        )
        state.log.append(
            {
                "cycle": state.cycles_done + 1,
                "ok": False,
                "events": events,
                "n_alta": state.n_alta,
                "n_baja": state.n_baja,
                "surplus": state.surplus,
            }
        )
        return False

    # conmutaciones: enganchar 'need' pesos (2 pernos cada uno ~ 2*E_perno o 1 ciclo make-before-break)
    for _ in range(need):
        state.n_alta -= 1
        state.n_on_obj += 1
        state.E_in += 2 * p.E_perno_J  # R y O
        events.append("enganche peso alta->objeto")

    if not state.sinks(p):
        events.append("FAIL: tras enganche aun no se hunde")
        state.log.append({"cycle": state.cycles_done + 1, "ok": False, "events": events})
        return False

    # --- 2) bajar ---
    # Masa total (para dinamica / log)
    m_down = p.m_obj_base + state.n_on_obj * p.m_peso
    W_pe_total = mgh(m_down, p.delta_h)
    # Contabilidad HONESTA de generacion electrica:
    # Solo el lastre EXTRA que NO vuelve con el objeto (n_sink - n_float pesos)
    # aporta PE "gastable". El resto se anula al flotar/subir el objeto.
    n_extra = max(0, state.n_on_obj - p.n_float)
    m_ballast = n_extra * p.m_peso
    W_pe_ballast = mgh(m_ballast, p.delta_h)
    W_drag_b = p.drag_frac * W_pe_ballast
    W_rec = p.eta_gen * (W_pe_ballast - W_drag_b)  # honest E_out este ciclo
    W_rec_naive = p.eta_gen * (W_pe_total * (1.0 - p.drag_frac))
    state.E_out += W_rec
    state.h_obj = "baja"
    events.append(
        f"baja m={m_down:.1f} kg  W_rec_honest={W_rec:.1f} J  "
        f"(naive_total={W_rec_naive:.1f} J)"
    )

    # --- 3) soltar extras en cota baja (dejar objeto en n_float para flotar) ---
    to_drop = state.n_on_obj - p.n_float
    if to_drop < 0:
        to_drop = 0
    for _ in range(to_drop):
        state.n_on_obj -= 1
        state.n_baja += 1
        state.E_in += 2 * p.E_perno_J
        events.append("perneo peso objeto->recorrido BAJA")

    if not state.floats(p):
        events.append("WARN: tras soltar no flota (n_on_obj=%d)" % state.n_on_obj)

    # --- 4) subir objeto por flotacion (coste electrico 0 en este modelo) ---
    state.h_obj = "alta"
    events.append("sube objeto ligero (boyancia) + recarga distancia en lado S")

    # --- 5) opcional: reponer 1 peso baja->alta (cierra almacen) ---
    if do_lift_reset and state.n_baja > 0:
        state.n_baja -= 1
        state.n_alta += 1
        W_lift = mgh(p.m_peso, p.delta_h) / p.eta_lift
        state.E_in += W_lift
        events.append(f"LIFT peso baja->alta  W_in={W_lift:.1f} J")

    state.cycles_done += 1
    state.log.append(
        {
            "cycle": state.cycles_done,
            "ok": True,
            "events": events,
            "n_alta": state.n_alta,
            "n_baja": state.n_baja,
            "n_on_obj": state.n_on_obj,
            "E_out": state.E_out,
            "E_in": state.E_in,
            "surplus": state.surplus,
            "W_rec_this": W_rec,
            "W_rec_naive_total": W_rec_naive,
            "m_ballast": m_ballast,
            "impar": state.cycles_done % 2 == 1,
        }
    )
    return True


def simulate(
    max_cycles: int = 12,
    p: Params | None = None,
    do_lift_reset: bool = False,
    label: str = "",
) -> State:
    p = p or Params()
    # reposo: n_float pesos en objeto (flota), stock en guia
    st = State(
        n_alta=p.n_stock_alta,
        n_baja=p.n_stock_baja,
        n_on_obj=p.n_float,
        h_obj="alta",
    )
    for _ in range(max_cycles):
        ok = run_cycle(st, p, do_lift_reset=do_lift_reset)
        if not ok:
            break
    st.log.append(
        {
            "summary": True,
            "label": label or ("con_lift" if do_lift_reset else "sin_lift"),
            "cycles_done": st.cycles_done,
            "n_alta": st.n_alta,
            "n_baja": st.n_baja,
            "E_out": st.E_out,
            "E_in": st.E_in,
            "surplus": st.surplus,
            "stopped_reason": "stock_alta_agotado_o_max"
            if st.cycles_done < max_cycles
            else "max_cycles",
        }
    )
    return st


def honest_surplus_from_log(st: State, p: Params) -> list[dict]:
    """Filas por ciclo desde el log (E_out ya es honest ballast-only)."""
    rows = []
    for e in st.log:
        if not e.get("ok"):
            continue
        w_rec = e["W_rec_this"]
        w_pernos = 4 * p.E_perno_J
        rows.append(
            {
                "cycle": e["cycle"],
                "W_rec_honest": w_rec,
                "W_rec_naive_total": e.get("W_rec_naive_total", 0.0),
                "W_pernos": w_pernos,
                "surplus_honest_cycle": w_rec - w_pernos,
                "n_alta": e["n_alta"],
                "n_baja": e["n_baja"],
                "impar": e["impar"],
            }
        )
    return rows


def main():
    p = Params()
    print("=" * 70)
    print("SIM: 3 flota / 4 se hunde | recarga SIEMPRE lado alto | stock finito")
    print("=" * 70)
    print(
        f"m_peso={p.m_peso} kg  Dh={p.delta_h} m  stock_alta={p.n_stock_alta}  "
        f"eta_gen={p.eta_gen}  m_g_h_peso={mgh(p.m_peso, p.delta_h):.1f} J"
    )
    print(f"n_float={p.n_float}  n_sink={p.n_sink}  n_puntas={p.n_puntas}")
    print()

    # A) Sin reponer pesos (almacen se agota) — el "impar XD"
    st_a = simulate(max_cycles=20, p=p, do_lift_reset=False, label="sin_lift_stock_finito")
    # B) Con reponer 1 peso por ciclo (cierra almacen)
    st_b = simulate(max_cycles=20, p=p, do_lift_reset=True, label="con_lift_cada_ciclo")

    def print_run(st: State, title: str):
        print("-" * 70)
        print(title)
        print(
            f"{'cic':>4} {'impar':>5} {'alta':>5} {'baja':>5} "
            f"{'E_out':>10} {'E_in':>10} {'surplus':>10} {'W_rec':>10}"
        )
        for e in st.log:
            if e.get("summary"):
                continue
            if not e.get("ok"):
                print(
                    f"{e.get('cycle', '?'):>4}  STOP  alta={e.get('n_alta')} "
                    f"baja={e.get('n_baja')}  | {e.get('events', [''])[0][:50]}"
                )
                continue
            print(
                f"{e['cycle']:4d} {'SI' if e['impar'] else 'NO':>5} "
                f"{e['n_alta']:5d} {e['n_baja']:5d} "
                f"{e['E_out']:10.1f} {e['E_in']:10.1f} {e['surplus']:10.1f} "
                f"{e['W_rec_this']:10.1f}"
            )
        sm = st.log[-1]
        print(
            f"\n  RESUMEN: ciclos={sm['cycles_done']}  alta={sm['n_alta']}  baja={sm['n_baja']}  "
            f"E_out={sm['E_out']:.1f}  E_in={sm['E_in']:.1f}  surplus_naive={sm['surplus']:.1f} J"
        )
        print(f"  razon fin: {sm.get('stopped_reason')}")

    print_run(st_a, "A) SIN subir pesos (stock alta se gasta) — efecto impar/inventario")
    print_run(st_b, "B) CON lift 1 peso/ciclo (cierra almacen) — peaje visible")

    # Contabilidad honesta en A
    honest = honest_surplus_from_log(st_a, p)
    print("\n" + "-" * 70)
    print("CONTABILIDAD HONESTA (solo lastre que se queda abajo; boyancia anula el resto)")
    print(f"{'cic':>4} {'impar':>5} {'rec_ballast':>12} {'pernos':>8} {'s_hon_ciclo':>12} {'alta':>5} {'baja':>5}")
    acc_h = 0.0
    for r in honest:
        acc_h += r["surplus_honest_cycle"]
        print(
            f"{r['cycle']:4d} {'SI' if r['impar'] else 'NO':>5} "
            f"{r['W_rec_honest']:12.1f} {r['W_pernos']:8.1f} "
            f"{r['surplus_honest_cycle']:12.1f} {r['n_alta']:5d} {r['n_baja']:5d}"
        )
    # Al final, el stock baja tiene N pesos; reponerlos TODOS costaria:
    n_left_baja = st_a.n_baja
    W_reset_all = n_left_baja * mgh(p.m_peso, p.delta_h) / p.eta_lift
    print(f"\n  surplus_honest acumulado (sin reset final): {acc_h:.1f} J")
    print(f"  coste de subir TODOS los pesos aparcados abajo ({n_left_baja}): {W_reset_all:.1f} J")
    print(f"  saldo si cierras el almacen al final: {acc_h - W_reset_all:.1f} J")
    print(
        "  -> El 'impar XD' deja pesos abajo; la energia 'ganada' es esa PE de lastre. "
        "Al reponer, se cancela (y eta la deja negativa)."
    )

    # Plots
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    def series(st):
        cyc, sout, sin, alta, baja, sur = [], [], [], [], [], []
        for e in st.log:
            if e.get("ok"):
                cyc.append(e["cycle"])
                sout.append(e["E_out"])
                sin.append(e["E_in"])
                alta.append(e["n_alta"])
                baja.append(e["n_baja"])
                sur.append(e["surplus"])
        return (
            np.array(cyc),
            np.array(sout),
            np.array(sin),
            np.array(alta),
            np.array(baja),
            np.array(sur),
        )

    ca, outa, ina, alta_a, baja_a, sura = series(st_a)
    cb, outb, inb, alta_b, baja_b, surb = series(st_b)

    ax = axes[0, 0]
    ax.plot(ca, outa, "g-o", label="E_out (naive)")
    ax.plot(ca, ina, "r-o", label="E_in (pernos)")
    ax.plot(ca, sura, "b--o", label="surplus naive")
    ax.set_title("A) Sin lift — energia naive")
    ax.set_xlabel("ciclo")
    ax.set_ylabel("J")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    for c in ca:
        if c % 2 == 1:
            ax.axvline(c, color="orange", alpha=0.15, lw=6)

    ax = axes[0, 1]
    ax.step(ca, alta_a, "b-o", where="mid", label="stock ALTA")
    ax.step(ca, baja_a, "m-o", where="mid", label="stock BAJA")
    ax.set_title("A) Inventario pesos (impar deja mas abajo)")
    ax.set_xlabel("ciclo")
    ax.set_ylabel("n pesos")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(cb, outb, "g-o", label="E_out")
    ax.plot(cb, inb, "r-o", label="E_in (pernos+lift)")
    ax.plot(cb, surb, "k--o", label="surplus")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title("B) Con lift cada ciclo — peaje visible")
    ax.set_xlabel("ciclo")
    ax.set_ylabel("J")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    # barras: surplus naive al parar A vs saldo con reset final
    labels = [
        "A surplus\nnaive al parar",
        "A honest\nsin reset",
        "A honest\n- reset all",
        "B surplus\ncon lift",
    ]
    vals = [
        st_a.surplus,
        acc_h,
        acc_h - W_reset_all,
        st_b.surplus,
    ]
    colors = ["#e67e22", "#3498db", "#27ae60" if vals[2] < 0 else "#e74c3c", "#9b59b6"]
    ax.bar(labels, vals, color=colors)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title("Cierre contable: naive vs honesto vs lift")
    ax.set_ylabel("J")
    ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle(
        "Kilometro lastre 3/4 + pernos | recarga siempre lado alto | stock finito",
        fontsize=12,
    )
    fig.tight_layout()
    png = OUT / "sim_enjambre_pesos_impar.png"
    fig.savefig(png, dpi=140, bbox_inches="tight")
    plt.close(fig)

    report = {
        "params": asdict(p),
        "run_A_sin_lift": {
            "cycles": st_a.cycles_done,
            "n_alta": st_a.n_alta,
            "n_baja": st_a.n_baja,
            "E_out": st_a.E_out,
            "E_in": st_a.E_in,
            "surplus_naive": st_a.surplus,
            "log": [e for e in st_a.log if not e.get("summary")],
        },
        "run_B_con_lift": {
            "cycles": st_b.cycles_done,
            "n_alta": st_b.n_alta,
            "n_baja": st_b.n_baja,
            "E_out": st_b.E_out,
            "E_in": st_b.E_in,
            "surplus": st_b.surplus,
            "log": [e for e in st_b.log if not e.get("summary")],
        },
        "honest_A": {
            "per_cycle": honest,
            "acc_honest_no_reset": acc_h,
            "reset_all_cost": W_reset_all,
            "saldo_si_cierras_almacen": acc_h - W_reset_all,
        },
        "verdict_es": (
            "Sin lift: el sistema hace tantos ciclos como pesos en stock alta; "
            "el surplus naive crece porque cobras m g h de lastre que dejas abajo. "
            "Impar/par solo marca el inventario. Al reponer todos los pesos, el saldo "
            "honesto cierra en negativo (eta). Con lift cada ciclo, surplus < 0 siempre."
        ),
    }
    jpath = OUT / "sim_enjambre_pesos_impar.json"
    jpath.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("\n" + "=" * 70)
    print(f"Grafica: {png}")
    print(f"JSON:    {jpath}")
    print(report["verdict_es"])


if __name__ == "__main__":
    main()
