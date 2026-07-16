#!/usr/bin/env python3
"""
AntiPC CLI v0.2 — calculadora de conocimiento desde CMD.

  antipc hash archivo.bin
  antipc hash --text "Hola"
  antipc mdc factor 10403
  antipc wave
  antipc mechanical 2 3
  antipc boot [--metrics out.json]
  antipc mk hash --text "Hola"
  antipc network demo [--hubs 4] [--duration 3]
  antipc version
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from hash_engine import get_backend, k3_hash_buffer, source_origin
from plugins.k3_plugin import K3HashPlugin
from plugins.mdc_factor_plugin import MdcFactorPlugin
from providers.mechanical_provider import MechanicalProvider
from providers.wave_provider import WaveProvider
from runtime.kernel import FlowKernel
from runtime.microkernel import AntiPCMicroKernel
from runtime.modes import WaveMode
from runtime.plugin import PluginContext
from runtime.reference import Reference, ReferenceRecord, ReferenceState


VERSION = "0.2.0-cmd"


def _kernel_with_plugins(*, wave: bool = False, wave_host: str = "8.8.8.8") -> FlowKernel:
    if wave:
        k = FlowKernel(mode=WaveMode(), wave_provider=WaveProvider(host=wave_host))
        k.network_permission = True
    else:
        k = FlowKernel()
    k.register_plugin(K3HashPlugin())
    k.register_plugin(MdcFactorPlugin())
    return k


def cmd_version(_: argparse.Namespace) -> int:
    print(f"AntiPC CLI {VERSION}")
    print(f"Hash backend : {get_backend()}")
    print(f"K3 origin    : {source_origin()}")
    return 0


def cmd_hash(args: argparse.Namespace) -> int:
    if args.text:
        data = args.text.encode("utf-8")
        label = "stdin-text"
    elif args.file:
        path = Path(args.file)
        if not path.is_file():
            print(f"ERROR: no existe {path}", file=sys.stderr)
            return 1
        data = path.read_bytes()
        label = path.name
    else:
        print("ERROR: usa --file o --text", file=sys.stderr)
        return 1

    kernel = _kernel_with_plugins(wave=args.wave, wave_host=args.wave_host)
    ref_in = kernel.ingest_raw(data, label=label)
    plugin = kernel.plugins["K3_HASH"]

    t0 = time.perf_counter()
    sig1 = kernel.submit_operation("K3_HASH", [ref_in])
    kernel.run_until_idle()
    t1 = time.perf_counter()

    t2 = time.perf_counter()
    sig2 = kernel.submit_operation("K3_HASH", [ref_in])
    kernel.run_until_idle()
    t3 = time.perf_counter()

    st = kernel.stats()
    digest = k3_hash_buffer(data)
    print(f"Input     : {label} ({len(data)} bytes)")
    print(f"K3 digest : {digest:08X}")
    print(f"Backend   : {get_backend()}")
    print(f"1st run   : {int((t1 - t0) * 1000)} ms  sig={sig1[:16]}...")
    print(f"2nd run   : {int((t3 - t2) * 1000)} ms  (reuse)")
    print(f"Knowledge : hits {st['knowledge_hits']}/{st['knowledge_queries']}")
    if st["knowledge_hits"] > 0:
        print("REUSE     : OK — segunda vez sin recalcular ALU")
    if args.wave and st.get("wave_inference") is not None:
        print(f"WaveMode  : inferencia={st['wave_inference']:.4f}  deferred={st['wave_deferred']}")
    return 0


def cmd_mdc_visual(args: argparse.Namespace) -> int:
    from mdc_lib.mdc_visual import animate_trains, launch_gui

    if args.n < 4:
        print("ERROR: n debe ser >= 4", file=sys.stderr)
        return 1
    if args.gui:
        launch_gui(args.n, proper_only=args.proper)
        return 0
    out = args.output or f"output/mdc_trains_{args.n}.gif"
    path = animate_trains(
        args.n,
        proper_only=args.proper,
        save_path=out,
        show=args.show,
        fps=args.fps,
    )
    if path:
        print(f"Animación : {path}")
    return 0


def cmd_mdc_analyze(args: argparse.Namespace) -> int:
    from mdc_lib.mdc_analyze import analyze, format_report

    if args.n < 4:
        print("ERROR: n debe ser >= 4", file=sys.stderr)
        return 1
    r = analyze(args.n, d_aux=args.d)
    proper = getattr(args, "proper", False)
    print(format_report(r, d_aux=args.d, proper_only=proper))
    if args.json:
        from mdc_lib.mdc_analyze import _is_proper_pair

        pairs = r.union_both
        if proper:
            pairs = [p for p in pairs if _is_proper_pair(p[0], p[1], r.n)]
        out = {
            "n": r.n,
            "union_both": [{"s": s, "t": t, "product": s * t} for s, t in pairs],
            "train_x": [
                {"x": c.x, "y": c.y, "s": c.s, "t": c.t, "k": c.k}
                for c in r.train_x.collisions
            ],
            "train_y": [
                {"y": c.y, "x": c.x, "s": c.s, "t": c.t, "k": c.k}
                for c in r.train_y.collisions
            ],
            "toy_factor": r.toy_factor,
            "aux_step": r.aux_step,
        }
        Path(args.json).write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"JSON      : {args.json}")
    if getattr(args, "gif", None):
        from mdc_lib.mdc_visual import animate_trains

        path = animate_trains(
            args.n,
            proper_only=proper,
            save_path=args.gif,
            show=False,
            fps=8,
        )
        if path:
            print(f"GIF       : {path}")
    return 0


def cmd_mdc_factor(args: argparse.Namespace) -> int:
    n = args.n
    kernel = _kernel_with_plugins()
    payload = n.to_bytes(8, "little")
    ref_in = Reference.create(signature=f"mdc-in-{n}", payload=payload)
    kernel.references[ref_in.ref_id] = ReferenceRecord(ref_in, ReferenceState.READY)

    plugin = kernel.plugins["MDC_FACTOR"]
    ctx = PluginContext([ref_in], [payload])
    if not plugin.validate(ctx):
        print(f"ERROR: N={n} fuera de rango toy (ver mdc_lib)", file=sys.stderr)
        return 1

    t0 = time.perf_counter()
    ref_out = plugin.execute(ctx)
    ms = int((time.perf_counter() - t0) * 1000)
    result = json.loads(ref_out.payload.decode() if ref_out.payload else "{}")
    print(f"N         : {n}")
    print(f"Factor    : {result.get('factor')}")
    print(f"Tiempo    : {ms} ms")
    return 0


def cmd_wave(args: argparse.Namespace) -> int:
    wp = WaveProvider(host=args.host)
    print(f"WaveProvider v0.1 — host {args.host}")
    for i in range(args.samples):
        s = wp.sample()
        print(
            f"  [{i + 1}] rtt={s.latency_us} us  jitter={s.jitter_us} us  "
            f"inferencia={s.inference:.4f}"
        )
        time.sleep(0.05)
    print(f"entropy byte: 0x{wp.as_entropy_byte():02X}")
    return 0


def cmd_network_demo(args: argparse.Namespace) -> int:
    from network.bd_pipeline import run_network_demo

    out = Path(args.out) if args.out else Path("output") / "network_bd_demo.json"
    m = run_network_demo(
        port=args.port,
        duration_s=args.duration,
        hubs=args.hubs,
        packets=args.packets,
        batch_size=args.batch_size,
        spawn=not args.no_hubs,
        out_path=out,
    )
    print()
    print("=" * 72)
    print("  AntiPC network demo — B (slot ring) → D (Grafcet)")
    print("=" * 72)
    print(f"  UDP in        : {m.packets_in}")
    print(f"  Grafcet out   : {m.packets_grafcet}")
    print(f"  Drops         : {m.drops}")
    print(f"  User copies   : {m.memcpy_user_copies}  (1× recvinto/slot, sin cola)")
    print(f"  Cache hits    : {m.cache_hits}")
    print(f"  Throughput    : {m.throughput_pps:.0f} pkt/s")
    print(f"  Grafcet step  : {m.grafcet_step}")
    print(f"  JSON          : {out}")
    print("=" * 72)
    return 0


def cmd_mechanical(args: argparse.Namespace) -> int:
    mp = MechanicalProvider()
    product = mp.multiply_increment(args.a, args.b)
    print(f"MechanicalProvider v0.1")
    print(f"  a*b (regla log) = {product}")
    for d in args.deltas:
        v = mp.step(d)
        print(f"  step({d}) -> log_state={v:.6f}")
    return 0


def cmd_boot(args: argparse.Namespace) -> int:
    mk = AntiPCMicroKernel()
    mk.boot()
    mk.run_boot_kop()
    mk.create_knowledge(b"AntiPC-CMD")
    metrics = Path(args.metrics) if args.metrics else None
    mk.shutdown(metrics_path=metrics)
    return 0


def _read_payload(args: argparse.Namespace) -> tuple[bytes, str] | None:
    if args.text:
        return args.text.encode("utf-8"), "stdin-text"
    if args.file:
        path = Path(args.file)
        if not path.is_file():
            print(f"ERROR: no existe {path}", file=sys.stderr)
            return None
        return path.read_bytes(), path.name
    return None


def cmd_mk_hash(args: argparse.Namespace) -> int:
    payload_info = _read_payload(args)
    if payload_info is None:
        print("ERROR: usa --file o --text", file=sys.stderr)
        return 1
    data, label = payload_info

    mk = AntiPCMicroKernel()
    mk.boot(verbose=not args.quiet)
    mk.run_boot_kop(verbose=not args.quiet)

    t0 = time.perf_counter()
    d1 = mk.hash_payload(data, verbose=not args.quiet)
    t1 = time.perf_counter()
    d2 = mk.hash_payload(data, verbose=not args.quiet)
    t2 = time.perf_counter()

    metrics = Path(args.metrics) if args.metrics else None
    mk.shutdown(verbose=not args.quiet, metrics_path=metrics)

    if d1 is None:
        return 1
    print(f"Input     : {label} ({len(data)} bytes)")
    print(f"K3 digest : {d1}")
    print(f"1st run   : {int((t1 - t0) * 1000)} ms")
    print(f"2nd run   : {int((t2 - t1) * 1000)} ms  (reuse)")
    if d1 == d2 and mk.metrics.kop_reused > 0:
        print("REUSE     : OK — HashKOP sin recalcular ALU")
    if metrics:
        print(f"Metrics   : {metrics}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="antipc",
        description="AntiPC — calculadora de conocimiento (CMD)",
    )
    parser.add_argument("--version", action="store_true", help="Mostrar version")
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("version", help="Backend y version").set_defaults(func=cmd_version)
    p_boot = sub.add_parser("boot", help="Microkernel demo corta")
    p_boot.add_argument(
        "--metrics",
        metavar="FILE",
        help="Exportar metricas JSON al cerrar",
    )
    p_boot.set_defaults(func=cmd_boot)

    p_mk = sub.add_parser("mk", help="Microkernel KOPs")
    p_mk_sub = p_mk.add_subparsers(dest="mk_cmd", required=True)
    p_mk_hash = p_mk_sub.add_parser("hash", help="HashKOP con reuse")
    p_mk_hash.add_argument("file", nargs="?", help="Archivo a hashear")
    p_mk_hash.add_argument("--text", "-t", help="Texto en lugar de archivo")
    p_mk_hash.add_argument("--metrics", metavar="FILE", help="Exportar metricas JSON")
    p_mk_hash.add_argument("--quiet", "-q", action="store_true", help="Solo resumen")
    p_mk_hash.set_defaults(func=cmd_mk_hash)

    p_hash = sub.add_parser("hash", help="K3 hash con reuse")
    p_hash.add_argument("file", nargs="?", help="Archivo a hashear")
    p_hash.add_argument("--text", "-t", help="Texto en lugar de archivo")
    p_hash.add_argument(
        "--wave",
        action="store_true",
        help="WaveMode: prioriza reuse segun latencia red (R003)",
    )
    p_hash.add_argument("--wave-host", default="8.8.8.8", help="Host RTT para --wave")
    p_hash.set_defaults(func=cmd_hash)

    p_mdc = sub.add_parser("mdc", help="MDC factorizacion toy")
    p_mdc_sub = p_mdc.add_subparsers(dest="mdc_cmd", required=True)
    p_fac = p_mdc_sub.add_parser("factor", help="Factorizar N")
    p_fac.add_argument("n", type=int)
    p_fac.set_defaults(func=cmd_mdc_factor)

    p_an = p_mdc_sub.add_parser(
        "analyze",
        help="Dos trenes X/Y: colisiones enteras y union",
    )
    p_an.add_argument("n", type=int, help="Entero a analizar")
    p_an.add_argument("-d", type=int, default=None, help="Denominador recta auxiliar")
    p_an.add_argument("--json", metavar="FILE", help="Exportar resultado JSON")
    p_an.add_argument(
        "--proper",
        action="store_true",
        help="Solo factores propios (excluir 1×n)",
    )
    p_an.add_argument("--gif", metavar="FILE", help="Generar animacion GIF de trenes")
    p_an.set_defaults(func=cmd_mdc_analyze)

    p_vis = p_mdc_sub.add_parser("visual", help="Animar dos trenes X/Y + colisiones")
    p_vis.add_argument("n", type=int)
    p_vis.add_argument("-o", "--output", help="Guardar GIF (default output/mdc_trains_N.gif)")
    p_vis.add_argument("--fps", type=int, default=8)
    p_vis.add_argument("--show", action="store_true", help="Mostrar ventana matplotlib")
    p_vis.add_argument("--gui", action="store_true", help="Ventana Tkinter estatica")
    p_vis.add_argument("--proper", action="store_true", help="Solo colisiones propias")
    p_vis.set_defaults(func=cmd_mdc_visual)

    p_wave = sub.add_parser("wave", help="Inferencia latencia red/WiFi")
    p_wave.add_argument("--host", default="8.8.8.8")
    p_wave.add_argument("--samples", type=int, default=5)
    p_wave.set_defaults(func=cmd_wave)

    p_net = sub.add_parser("network", help="Demos UDP / Grafcet")
    p_net_sub = p_net.add_subparsers(dest="network_cmd", required=True)
    p_demo = p_net_sub.add_parser(
        "demo",
        help="B slot-ring UDP → D Grafcet (existencia + reuse)",
    )
    p_demo.add_argument("--port", type=int, default=3333)
    p_demo.add_argument("--duration", type=float, default=3.0)
    p_demo.add_argument("--hubs", type=int, default=4)
    p_demo.add_argument("--packets", type=int, default=20_000)
    p_demo.add_argument("--batch-size", type=int, default=32)
    p_demo.add_argument("--out", metavar="FILE", help="Exportar metricas JSON")
    p_demo.add_argument(
        "--no-hubs",
        action="store_true",
        help="No lanzar hubs (esperar emision externa)",
    )
    p_demo.set_defaults(func=cmd_network_demo)

    p_mech = sub.add_parser("mechanical", help="Regla mecanica toy")
    p_mech.add_argument("a", type=float)
    p_mech.add_argument("b", type=float)
    p_mech.add_argument("deltas", nargs="*", type=float, default=[0.1, 0.5, 1.0])
    p_mech.set_defaults(func=cmd_mechanical)

    args = parser.parse_args(argv)
    if args.version and args.command is None:
        return cmd_version(args)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())