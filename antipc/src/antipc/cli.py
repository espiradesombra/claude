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
  antipc game demo [--players 128] [--shards 4]
  antipc version
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

_ANTIPC_SRC = os.path.dirname(os.path.abspath(__file__))
if _ANTIPC_SRC not in sys.path:
    sys.path.insert(0, _ANTIPC_SRC)

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


VERSION = "0.14.2-cmd"


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
    from native_engine import status_report

    print(f"AntiPC CLI {VERSION}")
    print(f"Hash backend : {get_backend()}")
    print(f"K3 origin    : {source_origin()}")
    print(status_report())
    return 0


def cmd_native_status(_: argparse.Namespace) -> int:
    from native_engine import status_report

    print(status_report())
    return 0


def cmd_native_bench(args: argparse.Namespace) -> int:
    from native_engine import bench_native

    r = bench_native(limit=args.limit, mdc_n=args.mdc_n)
    print("AntiPC native bench")
    print(f"  Backend       : {r['backend']}")
    print(f"  Sieve C       : {r['sieve_c_ms']} ms  count={r['sieve_count']}")
    print(f"  Sieve Python  : {r['sieve_py_ms']} ms  count={r['sieve_py_count']}")
    print(f"  K3 hash       : {r['hash_ms']} ms  {r['hash_hex']}")
    print(f"  MDC factor    : {r['mdc_ms']} ms  n={r['mdc_n']} -> {r['mdc_factor']}")
    if "geo_perimeter" in r:
        print(f"  Geo converge  : {r['geo_ms']} ms  perimeter={r['geo_perimeter']}")
    if "mdc_trains_ms" in r:
        print(f"  MDC trains C  : {r['mdc_trains_ms']} ms")
    if "sieve_hibrida_ms" in r:
        print(f"  Criba hibrida : {r['sieve_hibrida_ms']} ms  count={r.get('sieve_hibrida_count')}")
    if "newton_j" in r:
        print(f"  Newton(121)   : {r['newton_ms']} ms  j={r['newton_j']:.6f} ok={r.get('newton_ok')}")
    if "ksweep_factor" in r:
        print(f"  K-sweep       : {r['ksweep_ms']} ms  D={r['ksweep_factor']} evals={r.get('ksweep_evals')}")
    if "aleatorovix_ok" in r:
        print(f"  Aleatorovix   : {r['aleatorovix_ms']} ms  roundtrip={r['aleatorovix_ok']}")
    print(f"  Criva(10k)    : {r['criva_10k']:.6f}")
    if args.json:
        Path(args.json).write_text(json.dumps(r, indent=2), encoding="utf-8")
        print(f"  JSON          : {args.json}")
    return 0


def cmd_criba(args: argparse.Namespace) -> int:
    from native_engine import (
        get_backend,
        sieve_count,
        sieve_desmemoriada_count,
        sieve_hibrida_count,
        sieve_modular6k_count,
        sieve_primes,
    )

    t0 = time.perf_counter()
    if args.desmemoriada:
        count = sieve_desmemoriada_count(args.limit)
        label = "Criba desmemoriada VMA"
    elif args.modular6k:
        count = sieve_modular6k_count(args.limit)
        label = "Criba modular 6k±1"
    elif args.hibrida:
        count = sieve_hibrida_count(args.limit)
        label = "Criba híbrida VMA"
    else:
        count = sieve_count(args.limit)
        label = "Criba Eratóstenes"
    ms = int((time.perf_counter() - t0) * 1000)
    print(f"{label} ≤ {args.limit}")
    print(f"  Primos : {count}")
    print(f"  Tiempo : {ms} ms")
    print(f"  Backend: {get_backend()}")
    if args.list:
        primes = sieve_primes(args.limit, cap=min(count, args.max_list))
        print(f"  Lista  : {primes[: args.max_list]}")
    return 0


def cmd_newton(args: argparse.Namespace) -> int:
    from native_engine import newton_log

    t0 = time.perf_counter()
    r = newton_log(
        args.E,
        b=args.base,
        familia=args.familia,
        n_exp=args.n_exp,
        k_known=args.k,
    )
    ms = int((time.perf_counter() - t0) * 1000)
    print(f"Newton Rápido log_{args.base}({args.E})")
    print(f"  Familia   : {args.familia}")
    print(f"  j         : {r['j']:.12f}")
    print(f"  j exacto  : {r['j_exacto']:.12f}")
    print(f"  Error     : {r['error']:.2e}")
    print(f"  Iter      : {r['iteraciones']}")
    print(f"  Converged : {r['converged']}")
    print(f"  Backend   : {r.get('backend', 'python')}")
    print(f"  Tiempo    : {ms} ms")
    return 0


def cmd_mdc_ksweep(args: argparse.Namespace) -> int:
    from native_engine import get_backend, mdc_ksweep

    t0 = time.perf_counter()
    factor, evals = mdc_ksweep(args.n, predict=not args.classic)
    ms = int((time.perf_counter() - t0) * 1000)
    mode = "clásico" if args.classic else "predictivo"
    print(f"MDC K-sweep {mode} — n={args.n}")
    print(f"  Factor D  : {factor}")
    if evals:
        print(f"  Evals     : {evals}")
    print(f"  Tiempo    : {ms} ms")
    print(f"  Backend   : {get_backend()}")
    return 0


def cmd_mdc_jerk(args: argparse.Namespace) -> int:
    from mdc_lib.mdc_jerk import analyze_jerk, format_jerk_report

    if args.n < 4:
        print("ERROR: n debe ser >= 4", file=sys.stderr)
        return 1
    t0 = time.perf_counter()
    r = analyze_jerk(args.n, m=args.m)
    ms = int((time.perf_counter() - t0) * 1000)
    print(format_jerk_report(r))
    print(f"  Tiempo    : {ms} ms")
    if args.factorize:
        _run_mdc_v23_factor(args.n)
    return 0


def _run_mdc_v23_factor(n: int) -> None:
    """Factorización vía mdc_v23.py (filestot l5 / DeepSeek)."""
    import importlib.util

    l5 = Path(__file__).resolve().parents[3] / "filestot l5" / "mdc_v23.py"
    if not l5.is_file():
        print(f"  (mdc_v23 no encontrado: {l5})")
        return
    spec = importlib.util.spec_from_file_location("mdc_v23", l5)
    if spec is None or spec.loader is None:
        return
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if not hasattr(mod, "mdc_v23"):
        return
    print("\n--- mdc_v23 (DeepSeek L5) ---")
    factor, t_ms = mod.mdc_v23(n, verbose=True)
    if factor:
        print(f"  Factor    : {factor}  q={n // factor}")
    print(f"  Tiempo    : {t_ms:.2f} ms")


def cmd_mrauv_calibrar(args: argparse.Namespace) -> int:
    from mdc_lib.mrauv import calibrar, format_calibrar

    if args.n0 < 3:
        print("ERROR: n0 debe ser >= 3", file=sys.stderr)
        return 1
    t0 = time.perf_counter()
    r = calibrar(args.n0, dn=args.dn)
    ms = int((time.perf_counter() - t0) * 1000)
    print(format_calibrar(r))
    print(f"  Tiempo    : {ms} ms")
    return 0


def cmd_mrauv_densidad(args: argparse.Namespace) -> int:
    from mdc_lib.mrauv import densidad_report, format_densidad

    if args.n < 1:
        print("ERROR: n debe ser >= 1", file=sys.stderr)
        return 1
    t0 = time.perf_counter()
    r = densidad_report(args.n)
    ms = int((time.perf_counter() - t0) * 1000)
    print(format_densidad(r))
    print(f"  Tiempo    : {ms} ms")
    return 0


def cmd_discriminant_factor(args: argparse.Namespace) -> int:
    from mdc_lib.discriminant import discriminant_factor, format_result

    if args.n < 2:
        print("ERROR: n debe ser >= 2", file=sys.stderr)
        return 1
    t0 = time.perf_counter()
    r = discriminant_factor(args.n)
    ms = int((time.perf_counter() - t0) * 1000)
    print(format_result(r))
    print(f"  Tiempo    : {ms} ms")
    return 0


def cmd_discriminant_trajectory(args: argparse.Namespace) -> int:
    from mdc_lib.discriminant import delta_trajectory, format_trajectory

    if args.n < 2:
        print("ERROR: n debe ser >= 2", file=sys.stderr)
        return 1
    t0 = time.perf_counter()
    rows = delta_trajectory(args.n, show_all=args.all)
    ms = int((time.perf_counter() - t0) * 1000)
    print(format_trajectory(args.n, rows))
    print(f"  Tiempo    : {ms} ms")
    return 0


def cmd_mrauv_goldbach(args: argparse.Namespace) -> int:
    from mdc_lib.mrauv import (
        D,
        count_goldbach,
        format_goldbach_scan,
        scan_goldbach,
        sieve_primes,
    )

    if args.n_max < args.n_start:
        print("ERROR: --n-max debe ser >= --n-start", file=sys.stderr)
        return 1
    t0 = time.perf_counter()
    rows, all_ok = scan_goldbach(args.n_start, args.n_max, args.delta)
    print(format_goldbach_scan(rows, all_ok))

    if args.compare:
        max_two = max(2 * n for n in args.compare)
        primes = sieve_primes(max_two)
        pset = set(primes)
        print("\nComparación MRAUV vs Goldbach exacto:")
        print(f"{'2n':>10}  {'G_exact':>9}  {'D(n)·2n':>10}")
        print("-" * 35)
        for n in args.compare:
            two_n = 2 * n
            g = count_goldbach(two_n, pset)
            pred = D(n) * two_n
            print(f"{two_n:>10}  {g:>9}  {pred:>10.1f}")

    ms = int((time.perf_counter() - t0) * 1000)
    print(f"  Tiempo    : {ms} ms")
    return 0 if all_ok else 1


def cmd_geo_masivo(args: argparse.Namespace) -> int:
    from native_engine import geo_masivo_crypt, get_backend

    if args.text:
        data = args.text.encode("utf-8")
        label = "texto"
    elif args.file:
        path = Path(args.file)
        if not path.is_file():
            print(f"ERROR: no existe {path}", file=sys.stderr)
            return 1
        data = path.read_bytes()
        label = path.name
    else:
        data = args.text.encode("utf-8") if args.text else b""
        if not data:
            print("ERROR: usa --text o --file", file=sys.stderr)
            return 1

    semilla = args.semilla & 0xFFFF
    t0 = time.perf_counter()
    cifrado, hash_fase = geo_masivo_crypt(data, semilla)
    t1 = time.perf_counter()
    claro, _ = geo_masivo_crypt(cifrado, semilla, decrypt=True, hash_fase=hash_fase)
    t2 = time.perf_counter()

    ok = claro == data
    print("Aleatorovix + fase geométrica (C)")
    print(f"  Entrada    : {label} ({len(data)} bytes)")
    print(f"  Semilla    : {semilla}")
    print(f"  Hash fase  : {hash_fase:.12e}")
    print(f"  Cifrar     : {int((t1 - t0) * 1000)} ms")
    print(f"  Descifrar  : {int((t2 - t1) * 1000)} ms")
    print(f"  Round-trip : {'OK' if ok else 'FALLO'}")
    if ok and data.decode("utf-8", errors="replace").isprintable():
        print(f"  Texto      : {claro.decode('utf-8')}")
    print(f"  Backend    : {get_backend()}")
    return 0 if ok else 1


def cmd_geo(args: argparse.Namespace) -> int:
    from native_engine import geo_converge, get_backend, is_full_native

    if args.demo:
        bits = "0101101011000010"
        tales_s = "3,5,8,13,21"
        puntos_s = "6,12,18"
    else:
        bits = args.bits
        tales_s = args.tales
        puntos_s = args.puntos
    tales = [int(x) for x in tales_s.split(",")]
    puntos = [int(x) for x in puntos_s.split(",")]
    if not is_full_native():
        print("ERROR: requiere antipc_native.dll (scripts\\21_build_antipc_native.bat)", file=sys.stderr)
        return 1
    t0 = time.perf_counter()
    per = geo_converge(bits, tales, puntos)
    ms = int((time.perf_counter() - t0) * 1000)
    print(f"Convergencia geométrica (C)")
    print(f"  Bits      : {bits}")
    print(f"  Tales     : {tales}")
    print(f"  Puntos    : {puntos}")
    print(f"  Perímetro : {per / 1000:.3f}  (raw={per})")
    print(f"  Tiempo    : {ms} ms")
    print(f"  Backend   : {get_backend()}")
    return 0


def cmd_k3_hash(args: argparse.Namespace) -> int:
    payload_info = _read_payload(args)
    if payload_info is None:
        print("ERROR: usa --file o --text", file=sys.stderr)
        return 1
    data, label = payload_info
    t0 = time.perf_counter()
    digest = k3_hash_buffer(data, seed=args.seed)
    ms = int((time.perf_counter() - t0) * 1000)
    print(f"Input     : {label} ({len(data)} bytes)")
    print(f"K3 digest : {digest:08X}")
    print(f"Tiempo    : {ms} ms")
    print(f"Backend   : {get_backend()}")
    return 0


def cmd_k3_file(args: argparse.Namespace) -> int:
    from native_engine import k3_fingerprint_file, k3_hash_file

    path = Path(args.file)
    if not path.is_file():
        print(f"ERROR: no existe {path}", file=sys.stderr)
        return 1
    t0 = time.perf_counter()
    if args.fingerprint:
        size, digest = k3_fingerprint_file(path)
    else:
        size = path.stat().st_size
        digest = k3_hash_file(path, seed=args.seed)
    ms = int((time.perf_counter() - t0) * 1000)
    print(f"Archivo   : {path.name} ({size} bytes)")
    print(f"K3 digest : {digest:08X}")
    print(f"Tiempo    : {ms} ms (streaming C)")
    print(f"Backend   : {get_backend()}")
    return 0


def cmd_k3_verify(_: argparse.Namespace) -> int:
    from verify_k3 import main as verify_main

    return int(verify_main())


def cmd_k3_dedup(args: argparse.Namespace) -> int:
    from native_engine import k3_dedup_paths

    paths: list[str] = list(args.files or [])
    if args.dir:
        root = Path(args.dir)
        if not root.is_dir():
            print(f"ERROR: no existe directorio {root}", file=sys.stderr)
            return 1
        paths.extend(str(p) for p in root.rglob("*") if p.is_file())
    if not paths:
        print("ERROR: indica archivos o --dir", file=sys.stderr)
        return 1

    t0 = time.perf_counter()
    report = k3_dedup_paths(paths)
    ms = int((time.perf_counter() - t0) * 1000)

    print(f"k3dedup — {report['files_analyzed']} archivos, {report['duplicate_groups']} grupos")
    for i, g in enumerate(report["groups"], 1):
        print(f"=== Grupo {i} (tam={g['size']} bytes, hash=0x{g['hash']}, {g['copies']} copias) ===")
        for p in g["paths"]:
            print(f"  {p}")
    print(f"~{report['bytes_recoverable']} bytes recuperables")
    print(f"Tiempo    : {ms} ms")
    print(f"Backend   : {get_backend()}")
    if args.json:
        Path(args.json).write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"JSON      : {args.json}")
    return 0


def cmd_k3_redundant(args: argparse.Namespace) -> int:
    from native_engine import k3_redundant_hashes

    payload_info = _read_payload(args)
    if payload_info is None:
        print("ERROR: usa --file o --text", file=sys.stderr)
        return 1
    data, label = payload_info
    hashes = k3_redundant_hashes(data, replicas=args.replicas)
    print(f"Input     : {label} ({len(data)} bytes)")
    for i, h in enumerate(hashes):
        print(f"  canal[{i}] : {h:08X}")
    print(f"Backend   : {get_backend()}")
    return 0


def cmd_k3_heavy(args: argparse.Namespace) -> int:
    from native_engine import k3_heavy_hash

    payload_info = _read_payload(args)
    if payload_info is None:
        print("ERROR: usa --file o --text", file=sys.stderr)
        return 1
    data, label = payload_info
    t0 = time.perf_counter()
    digest = k3_heavy_hash(data)
    ms = int((time.perf_counter() - t0) * 1000)
    print(f"Input     : {label} ({len(data)} bytes)")
    print(f"Heavy K3  : {digest:08X}")
    print(f"Tiempo    : {ms} ms")
    print(f"Backend   : {get_backend()}")
    return 0


def cmd_k3_hamming(args: argparse.Namespace) -> int:
    from native_engine import k3_hamming, k3_similarity

    a = int(args.hash_a, 16)
    b = int(args.hash_b, 16)
    dist = k3_hamming(a, b)
    sim = k3_similarity(a, b)
    print(f"Hash A    : {a:08X}")
    print(f"Hash B    : {b:08X}")
    print(f"Hamming   : {dist} bits")
    print(f"Similitud : {sim:.4f}")
    print(f"Backend   : {get_backend()}")
    return 0


def cmd_k3_simil(args: argparse.Namespace) -> int:
    from k3_text_tools import k3_simil_pairs

    paths = list(args.files or [])
    if args.dir:
        root = Path(args.dir)
        paths.extend(str(p) for p in root.rglob("*.txt") if p.is_file())
    if len(paths) < 2:
        print("ERROR: indica >=2 archivos o --dir con .txt", file=sys.stderr)
        return 1
    pairs = k3_simil_pairs(paths, threshold=args.threshold)
    print(f"k3simil — {len(paths)} docs, umbral {args.threshold:.0%}")
    for p in pairs:
        print(f"  {p['similarity']*100:5.1f}%  {p['path_a']}")
        print(f"         <->  {p['path_b']}")
    if not pairs:
        print("  (ningún par superó el umbral)")
    return 0


def cmd_k3_search_index(args: argparse.Namespace) -> int:
    from k3_text_tools import k3_search_index

    paths: list[str] = list(args.files or [])
    if args.dir:
        root = Path(args.dir)
        for ext in ("*.txt", "*.md", "*.py"):
            paths.extend(str(p) for p in root.rglob(ext) if p.is_file())
    if not paths:
        print("ERROR: indica archivos o --dir", file=sys.stderr)
        return 1
    r = k3_search_index(paths, args.output)
    print(f"k3search index — {r['documents']} docs, {r['terms']} términos")
    print(f"  Índice: {r['index']}")
    return 0


def cmd_k3_search_query(args: argparse.Namespace) -> int:
    from k3_text_tools import k3_search_query

    hits = k3_search_query(args.index, args.words, top=args.top)
    print(f"k3search query — {len(hits)} resultados")
    for h in hits:
        print(f"  {h['score']:6.1f}  {h['path']}")
    return 0


def cmd_k3_stream_xor(args: argparse.Namespace) -> int:
    from native_engine import k3_stream_xor

    payload_info = _read_payload(args)
    if payload_info is None:
        print("ERROR: usa --file o --text", file=sys.stderr)
        return 1
    data, label = payload_info
    buf = bytearray(data)
    k3_stream_xor(buf, base=args.base, rel=args.rel)
    print(f"Input     : {label} ({len(data)} bytes)")
    print(f"XOR K3    : base={args.base} rel={args.rel}")
    print(f"Hex out   : {buf.hex()}")
    print(f"Backend   : {get_backend()}")
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


def cmd_game_demo(args: argparse.Namespace) -> int:
    from network.game_server import run_game_demo

    out = Path(args.out) if args.out else Path("output") / "game_mmo_demo.json"
    m, world = run_game_demo(
        port=args.port,
        duration_s=args.duration,
        players=args.players,
        tick_hz=args.tick_hz,
        batch_size=args.batch_size,
        num_shards=args.shards,
        duplicate_ratio=args.duplicate_ratio,
        spawn_bots_flag=not args.no_bots,
        out_path=out,
    )
    print()
    print("=" * 72)
    print("  AntiPC MMO — B (slot-ring) → D (Grafcet) → WorldState")
    print("=" * 72)
    print(f"  Puerto        : {m.port}")
    print(f"  Jugadores     : {m.players}")
    print(f"  UDP in        : {m.packets_in}")
    print(f"  Movimientos   : {m.moves_applied}  ({m.moves_per_sec:.0f}/s)")
    print(f"  Entidades     : {m.entities_active}")
    print(f"  Drops         : {m.drops}")
    print(f"  Existencia OK : {m.moves_applied}  rechazados={m.rejected_existence}")
    print(f"  Grafcet cache : {m.grafcet_cache_hits}")
    print(f"  Estado cache  : {m.state_cache_hits}  (paquetes duplicados)")
    print(f"  Throughput    : {m.throughput_pps:.0f} pkt/s")
    print(f"  Max tick      : {m.max_tick}")
    print(f"  Shards        : {m.num_shards}")
    if world:
        print("  Muestra mundo :")
        for row in world[:5]:
            print(f"    entity {row['entity']} shard={row['shard']} pos={row['pos']}")
    print(f"  JSON          : {out}")
    print("=" * 72)
    return 0


def cmd_game_cluster(args: argparse.Namespace) -> int:
    from network.game_server import run_game_cluster_demo

    out = Path(args.out) if args.out else Path("output") / "game_cluster_demo.json"
    hub_hosts = [h.strip() for h in args.hub_hosts.split(",") if h.strip()] if args.hub_hosts else None
    m, world = run_game_cluster_demo(
        port=args.port,
        duration_s=args.duration,
        players=args.players,
        tick_hz=args.tick_hz,
        batch_size=args.batch_size,
        num_shards=args.shards,
        hubs=args.hubs,
        hub_hosts=hub_hosts,
        hub_base_port=args.hub_base_port,
        duplicate_ratio=args.duplicate_ratio,
        spawn_bots_flag=not args.no_bots,
        out_path=out,
    )
    print()
    print("=" * 72)
    print("  AntiPC MMO CLUSTER — B → L3 hubs (WORK/RESULT) → D → WorldState")
    print("=" * 72)
    print(f"  Puerto        : {m.port}")
    print(f"  Jugadores     : {m.players}")
    print(f"  Hubs L3       : {m.hubs_authenticated}/{m.hubs} autenticados (HMAC)")
    print(f"  UDP in        : {m.packets_in}")
    print(f"  Movimientos   : {m.moves_applied}  ({m.moves_per_sec:.0f}/s)")
    print(f"  Entidades     : {m.entities_active}")
    print(f"  Remote heavy  : {m.remote_validations} enviados, {m.remote_validations_ok} OK")
    if m.remote_fallback:
        print(f"  Fallback local: {m.remote_fallback}")
    print(f"  Grafcet cache : {m.grafcet_cache_hits}")
    print(f"  Estado cache  : {m.state_cache_hits}")
    print(f"  Throughput    : {m.throughput_pps:.0f} pkt/s")
    print(f"  Shards        : {m.num_shards}")
    if world:
        print("  Muestra mundo :")
        for row in world[:5]:
            print(f"    entity {row['entity']} shard={row['shard']} pos={row['pos']}")
    print(f"  JSON          : {out}")
    print("=" * 72)
    return 0


def cmd_game_bots(args: argparse.Namespace) -> int:
    from network.game_bots import run_bot_swarm

    run_bot_swarm(
        master=args.master,
        port=args.port,
        players=args.players,
        duration_s=args.duration,
        tick_hz=args.tick_hz,
        duplicate_ratio=args.duplicate_ratio,
    )
    return 0


def cmd_libro_list(_: argparse.Namespace) -> int:
    from libros.bridge import format_libro_info, list_libros

    print()
    print("=" * 72)
    print("  Libros VMA 1–6 — métodos matemáticos → AntiPC")
    print("=" * 72)
    for info in list_libros():
        print()
        print(format_libro_info(info))
    print("=" * 72)
    return 0


def cmd_libro_info(args: argparse.Namespace) -> int:
    from libros.bridge import LIBRO_REGISTRY, format_libro_info

    n = int(args.numero)
    if n not in LIBRO_REGISTRY:
        print(f"Libro {n} no registrado (use 1–6)")
        return 1
    print(format_libro_info(LIBRO_REGISTRY[n]))
    return 0


def cmd_libro_metodos(args: argparse.Namespace) -> int:
    from libros.bridge import format_metodos_table

    libro = int(args.libro) if args.libro else None
    print(format_metodos_table(libro=libro))
    return 0


def cmd_libro_run(args: argparse.Namespace) -> int:
    from libros.bridge import format_libro_info, resolve_metodo, run_metodo

    libro = int(args.numero)
    try:
        info, metodo = resolve_metodo(libro, args.metodo)
    except (KeyError, FileNotFoundError) as exc:
        print(f"Error: {exc}")
        return 1

    print(format_libro_info(info))
    print()
    print(f"Ejecutando método: {metodo.id} — {metodo.nombre}")
    try:
        result = run_metodo(libro, metodo.id, extra_n=args.n)
    except RuntimeError as exc:
        print(f"No ejecutable: {exc}")
        return 1

    if result.stdout.strip():
        print(result.stdout)
    if result.stderr.strip():
        print(result.stderr, file=sys.stderr)
    print(f"Exit {result.returncode} · {result.elapsed_s:.1f}s")
    return 0 if result.returncode == 0 else 1


def cmd_gemelo_list(_: argparse.Namespace) -> int:
    from gemelos.bridge import format_gemelo_info, list_gemelos

    print()
    print("=" * 72)
    print("  AntiPC gemelos — ZypyZape · Quijote · Kilòmetre (repo VMA)")
    print("=" * 72)
    for info in list_gemelos():
        print()
        print(format_gemelo_info(info))
    print("=" * 72)
    return 0


def cmd_gemelo_info(args: argparse.Namespace) -> int:
    from gemelos.bridge import GEMELO_REGISTRY, format_gemelo_info

    key = args.name.lower()
    if key not in GEMELO_REGISTRY:
        print(f"Gemelo desconocido: {args.name}")
        print(f"Disponibles: {', '.join(GEMELO_REGISTRY)}")
        return 1
    print(format_gemelo_info(GEMELO_REGISTRY[key]))
    return 0


def cmd_gemelo_run(args: argparse.Namespace) -> int:
    from gemelos.bridge import format_gemelo_info, resolve_script, run_gemelo

    key = args.name.lower()
    try:
        info, script = resolve_script(key, args.variant)
    except (KeyError, FileNotFoundError) as exc:
        print(f"Error: {exc}")
        return 1

    extra: list[str] = list(args.extra or [])
    if key == "zypyzape" and script.variant == "viability":
        extra.extend(
            [
                "--hubs", str(args.hubs),
                "--packets", str(args.packets),
                "--duration", str(args.duration),
            ]
        )
        if args.out:
            extra.extend(["--out", args.out])
        if args.no_hubs:
            extra.append("--no-hubs")

    print(format_gemelo_info(info))
    print()
    print(f"Ejecutando: {script.path.name} [{script.variant}]")
    try:
        result = run_gemelo(key, variant=script.variant, extra_args=extra or None)
    except FileNotFoundError as exc:
        print(f"Error: {exc}")
        return 1
    except subprocess.TimeoutExpired:
        print("Error: tiempo de espera agotado")
        return 1

    if result.stdout.strip():
        print(result.stdout)
    if result.stderr.strip():
        print(result.stderr, file=sys.stderr)
    print(f"Exit {result.returncode} · {result.elapsed_s:.1f}s")
    return 0 if result.returncode == 0 else 1


def cmd_deepseek_list(_: argparse.Namespace) -> int:
    from deepseek.bridge import format_deepseek_info, format_deepseek_table, list_deepseek

    print(format_deepseek_table())
    print()
    for entry in list_deepseek():
        print(format_deepseek_info(entry))
        print()
    return 0


def cmd_deepseek_info(args: argparse.Namespace) -> int:
    from deepseek.bridge import DEEPSEEK_REGISTRY, format_deepseek_info

    key = args.key.lower()
    if key not in DEEPSEEK_REGISTRY:
        print(f"Entrada desconocida: {args.key}")
        print(f"Disponibles: {', '.join(DEEPSEEK_REGISTRY)}")
        return 1
    print(format_deepseek_info(DEEPSEEK_REGISTRY[key]))
    return 0


def cmd_deepseek_run(args: argparse.Namespace) -> int:
    from deepseek.bridge import format_deepseek_info, resolve_script, run_script

    key = args.key.lower()
    try:
        entry, script = resolve_script(key, args.script)
    except (KeyError, FileNotFoundError) as exc:
        print(f"Error: {exc}")
        return 1

    extra: list[str] = list(args.extra or [])
    if args.n is not None and key == "mdc-u":
        extra = ["-c", f"from mdc_v23 import mdc_v23; mdc_v23({args.n}, verbose=True)"]

    print(format_deepseek_info(entry))
    print()
    print(f"Ejecutando: {script.path.name} [{script.id}]")
    try:
        if extra and extra[0] == "-c":
            proc = subprocess.run(
                [sys.executable, "-c", extra[1]],
                cwd=str(script.path.parent),
                capture_output=True,
                text=True,
                timeout=180.0,
            )
            result_stdout = proc.stdout
            result_stderr = proc.stderr
            rc = proc.returncode
            elapsed = 0.0
        else:
            run_r = run_script(key, script.id, extra_args=extra or None)
            result_stdout = run_r.stdout
            result_stderr = run_r.stderr
            rc = run_r.returncode
            elapsed = run_r.elapsed_s
    except subprocess.TimeoutExpired:
        print("Error: tiempo de espera agotado")
        return 1

    if result_stdout.strip():
        print(result_stdout)
    if result_stderr.strip():
        print(result_stderr, file=sys.stderr)
    if elapsed:
        print(f"Exit {rc} · {elapsed:.1f}s")
    else:
        print(f"Exit {rc}")
    return 0 if rc == 0 else 1


def cmd_teorema_list(_: argparse.Namespace) -> int:
    from teoremas.bridge import format_teoremas_table

    print(format_teoremas_table())
    return 0


def cmd_teorema_info(args: argparse.Namespace) -> int:
    from teoremas.bridge import format_teorema_info, resolve_teorema

    try:
        if args.numero is not None:
            info = resolve_teorema(numero=int(args.numero))
        else:
            info = resolve_teorema(teorema_id=args.id)
    except (KeyError, ValueError) as exc:
        print(f"Error: {exc}")
        return 1
    print(format_teorema_info(info, body=not args.no_body))
    return 0


def cmd_network_bench(args: argparse.Namespace) -> int:
    from architecture.antipc_stack import AntiPCStack
    from network.network_compute import run_network_bench

    stack = AntiPCStack()
    if not stack.request_network_permission():
        print("Error: permiso de red HMAC denegado")
        return 1

    try:
        result = run_network_bench(
            stack.fabric,
            count=args.count,
            hubs=args.hubs,
            window=args.window,
            batch_size=args.batch,
            payload_size=args.payload,
            timeout_s=args.timeout,
        )
    finally:
        stack.stop_network()

    out = Path(args.out) if args.out else Path("output") / "network_bench.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "version": VERSION,
        "items": result.items,
        "hubs": result.hubs,
        "parity_ok": result.parity_ok,
        "speedup": round(result.speedup, 3),
        "classic": {
            "received": result.classic.received,
            "elapsed_s": round(result.classic.elapsed_s, 4),
            "throughput": round(result.classic.throughput, 1),
            "singles_sent": result.classic.singles_sent,
        },
        "pipelined": {
            "received": result.pipelined.received,
            "elapsed_s": round(result.pipelined.elapsed_s, 4),
            "throughput": round(result.pipelined.throughput, 1),
            "batches_sent": result.pipelined.batches_sent,
            "batches_recv": result.pipelined.batches_recv,
            "singles_sent": result.pipelined.singles_sent,
        },
    }
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print()
    print("=" * 72)
    print("  AntiPC network bench — clásico vs pipelined (WORK_BATCH)")
    print("=" * 72)
    print(f"  Items         : {result.items}")
    print(f"  Hubs          : {result.hubs}")
    print(f"  Paridad K3    : {'OK' if result.parity_ok else 'FALLO'}")
    print(f"  Clásico       : {result.classic.received} en {result.classic.elapsed_s:.3f}s "
          f"({result.classic.throughput:.0f} hash/s)")
    print(f"  Pipelined     : {result.pipelined.received} en {result.pipelined.elapsed_s:.3f}s "
          f"({result.pipelined.throughput:.0f} hash/s)")
    print(f"  Speedup       : {result.speedup:.2f}x")
    print(f"  Lotes UDP     : send={result.pipelined.batches_sent} recv={result.pipelined.batches_recv}")
    print(f"  JSON          : {out}")
    print("=" * 72)
    return 0 if result.parity_ok and result.pipelined.received == result.items else 1


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


def cmd_industrial_inventory(args: argparse.Namespace) -> int:
    from industrial.inventory import export_inventory, format_inventory_text

    out_txt = Path(args.out) if args.out else Path("output") / "inventario_industrial.txt"
    out_json = Path(args.json) if args.json else Path("output") / "inventario_industrial.json"
    inv = export_inventory(cli_version=VERSION, txt_path=out_txt, json_path=out_json)
    if not args.quiet:
        print(format_inventory_text(inv))
    print(f"TXT  : {out_txt}")
    print(f"JSON : {out_json}")
    ok = inv.dll_present and inv.health.get("k3_parity") == "OK"
    print(f"Estado industrial: {'OPERATIVO' if ok else 'REVISAR'}")
    return 0 if ok else 1


def cmd_industrial_audit(args: argparse.Namespace) -> int:
    """Demo industrial rápida + inventario."""
    from industrial.inventory import export_inventory, format_inventory_text
    from industrial.runtime import IndustrialRuntime
    from plugins.k3_dedup_plugin import K3DedupPlugin
    from plugins.k3_file_plugin import K3FilePlugin
    from runtime.profile import ExecutionProfile

    root = Path(args.dir) if args.dir else Path(__file__).resolve().parent
    files = [str(p) for p in root.rglob("*") if p.is_file()][: args.max_files]

    rt = IndustrialRuntime(profile=ExecutionProfile.industrial_audit())
    rt.bootstrap(K3FilePlugin(), K3DedupPlugin())
    t0 = time.perf_counter()
    refs = rt.scan_files(files, workers=args.workers)
    dedup = rt.run_dedup_report(refs)
    elapsed = int((time.perf_counter() - t0) * 1000)
    stats = rt.report()

    out_dir = Path(args.out_dir) if args.out_dir else Path("output")
    rt.telemetry.export_csv(out_dir / "telemetria_industrial_audit.csv")
    rt.telemetry.export_json(out_dir / "informe_industrial_audit.json", extra=stats)

    inv = export_inventory(
        cli_version=VERSION,
        txt_path=out_dir / "inventario_industrial.txt",
        json_path=out_dir / "inventario_industrial.json",
    )

    print("=" * 72)
    print("  AntiPC AUDITORÍA INDUSTRIAL")
    print("=" * 72)
    print(f"  Directorio   : {root}")
    print(f"  Ficheros     : {len(refs)}")
    print(f"  Tiempo       : {elapsed} ms")
    print(f"  ALU ejecutado: {stats.get('executed', 0)}")
    print(f"  Knowledge hit: {stats.get('knowledge_hits', 0)}")
    if dedup:
        print(f"  Dup grupos   : {dedup.metadata.get('duplicate_groups', 0)}")
        print(f"  Recuperable  : {dedup.metadata.get('bytes_recoverable', 0)} B")
    print(f"  Telemetría   : {out_dir}")
    print(f"  Estado       : {'OPERATIVO' if inv.dll_present else 'SIN DLL'}")
    if args.full:
        print()
        print(format_inventory_text(inv))
    print("=" * 72)
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

    p_ks = p_mdc_sub.add_parser("ksweep", help="K-sweep entero (C predictivo)")
    p_ks.add_argument("n", type=int)
    p_ks.add_argument("--classic", action="store_true", help="K-sweep clásico sin predicción")
    p_ks.set_defaults(func=cmd_mdc_ksweep)

    p_jk = p_mdc_sub.add_parser("jerk", help="MDC-U pinça 4+4 + invariante jerk")
    p_jk.add_argument("n", type=int, help="Entero a analizar")
    p_jk.add_argument("-m", type=int, default=None, help="Índice m (default ~√N)")
    p_jk.add_argument(
        "--factorize",
        action="store_true",
        help="Tras jerk, factorizar con mdc_v23 (DeepSeek L5)",
    )
    p_jk.set_defaults(func=cmd_mdc_jerk)

    p_mra = sub.add_parser("mrauv", help="MRAUV densidad cinemática (Libro 2)")
    p_mra_sub = p_mra.add_subparsers(dest="mrauv_cmd", required=True)

    p_mrc = p_mra_sub.add_parser("calibrar", help="Calibración 3 puntos V0,a0,j")
    p_mrc.add_argument("n0", type=int, help="Punto inicial n0")
    p_mrc.add_argument("--dn", type=int, default=None, help="Paso Δ (default 2√n0)")
    p_mrc.set_defaults(func=cmd_mrauv_calibrar)

    p_mrd = p_mra_sub.add_parser("densidad", help="L(n), m(n), D(n), criterio 2-primos")
    p_mrd.add_argument("n", type=int)
    p_mrd.set_defaults(func=cmd_mrauv_densidad)

    p_mrg = p_mra_sub.add_parser("goldbach", help="Barrido criterio MRAUV-Goldbach [HEUR]")
    p_mrg.add_argument("--n-start", type=int, default=1000)
    p_mrg.add_argument("--n-max", type=int, default=50_000)
    p_mrg.add_argument("--delta", type=int, default=5000)
    p_mrg.add_argument(
        "--compare",
        type=int,
        nargs="*",
        metavar="N",
        help="Comparar con conteo exacto en estos n (ej. 100 500 1000)",
    )
    p_mrg.set_defaults(func=cmd_mrauv_goldbach)

    p_dis = sub.add_parser("discriminant", help="Método discriminante Δ(S)=k² (Libro 5)")
    p_dis_sub = p_dis.add_subparsers(dest="discriminant_cmd", required=True)

    p_disf = p_dis_sub.add_parser("factor", help="Factorizar / filtro primo por cuadrado perfecto")
    p_disf.add_argument("n", type=int, help="Entero N a analizar")
    p_disf.set_defaults(func=cmd_discriminant_factor)

    p_dist = p_dis_sub.add_parser("trajectory", help="Trayectoria Δ(S) hasta parada determinista")
    p_dist.add_argument("n", type=int)
    p_dist.add_argument("--all", action="store_true", help="Mostrar todos los pasos S")
    p_dist.set_defaults(func=cmd_discriminant_trajectory)

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
    p_bench = p_net_sub.add_parser(
        "bench",
        help="Benchmark clásico vs pipelined (hubs WORK/RESULT + lotes UDP)",
    )
    p_bench.add_argument("--count", type=int, default=512, help="WORK items")
    p_bench.add_argument("--hubs", type=int, default=4)
    p_bench.add_argument("--window", type=int, default=128, help="Ventana en vuelo")
    p_bench.add_argument("--batch", type=int, default=24, help="Tamaño lote WORK_BATCH")
    p_bench.add_argument("--payload", type=int, default=64, help="Bytes por payload")
    p_bench.add_argument("--timeout", type=float, default=15.0)
    p_bench.add_argument("--out", metavar="FILE", help="Exportar metricas JSON")
    p_bench.set_defaults(func=cmd_network_bench)

    p_gem = sub.add_parser(
        "gemelo",
        help="Gemelos digitales — ZypyZape, Quijote, Kilòmetre (repo VMA)",
    )
    p_gem_sub = p_gem.add_subparsers(dest="gemelo_cmd", required=True)
    p_gem_sub.add_parser("list", help="Listar gemelos y rutas en repo").set_defaults(
        func=cmd_gemelo_list
    )
    p_gi = p_gem_sub.add_parser("info", help="Detalle matemático de un gemelo")
    p_gi.add_argument("name", choices=["zypyzape", "quijote", "kilometre"])
    p_gi.set_defaults(func=cmd_gemelo_info)
    p_gr = p_gem_sub.add_parser("run", help="Ejecutar script del gemelo (numpy/mpl)")
    p_gr.add_argument("name", choices=["zypyzape", "quijote", "kilometre"])
    p_gr.add_argument(
        "--variant", "-v",
        help="viability|emitter (zypyzape), v48|v5 (quijote), v15|v14 (kilometre)",
    )
    p_gr.add_argument("--hubs", type=int, default=4, help="Solo zypyzape viability")
    p_gr.add_argument("--packets", type=int, default=15_000)
    p_gr.add_argument("--duration", type=float, default=2.5)
    p_gr.add_argument("--out", metavar="FILE", help="JSON salida (zypyzape)")
    p_gr.add_argument("--no-hubs", action="store_true")
    p_gr.add_argument("extra", nargs="*", help="Args extra al script Python")
    p_gr.set_defaults(func=cmd_gemelo_run)

    p_lib = sub.add_parser("libro", help="Métodos Libros 1–6 VMA → AntiPC")
    p_lib_sub = p_lib.add_subparsers(dest="libro_cmd", required=True)
    p_lib_sub.add_parser("list", help="Los 6 libros y sus métodos").set_defaults(
        func=cmd_libro_list
    )
    p_li = p_lib_sub.add_parser("info", help="Detalle de un libro (1–6)")
    p_li.add_argument("numero", type=int, choices=[1, 2, 3, 4, 5, 6])
    p_li.set_defaults(func=cmd_libro_info)
    p_lm = p_lib_sub.add_parser("metodos", help="Tabla método → comando AntiPC")
    p_lm.add_argument("--libro", type=int, choices=[1, 2, 3, 4, 5, 6], help="Filtrar libro")
    p_lm.set_defaults(func=cmd_libro_metodos)
    p_lr = p_lib_sub.add_parser("run", help="Ejecutar método con CLI AntiPC")
    p_lr.add_argument("numero", type=int, choices=[1, 2, 3, 4, 5, 6])
    p_lr.add_argument("--metodo", "-m", help="ID método (ej. L5-ksweep)")
    p_lr.add_argument("--n", type=int, help="Sustituir N en mdc factor/analyze/ksweep")
    p_lr.set_defaults(func=cmd_libro_run)

    p_ds = sub.add_parser("deepseek", help="Corpus DeepSeek 6 2026 (MDC-U, L5, tablas)")
    p_ds_sub = p_ds.add_subparsers(dest="deepseek_cmd", required=True)
    p_ds_sub.add_parser("list", help="Tabla + entradas registradas").set_defaults(
        func=cmd_deepseek_list
    )
    p_dsi = p_ds_sub.add_parser("info", help="Detalle de una entrada")
    p_dsi.add_argument("key", choices=["mdc-u", "libros-tabla", "criba-aleatorovix"])
    p_dsi.set_defaults(func=cmd_deepseek_info)
    p_dsr = p_ds_sub.add_parser("run", help="Ejecutar script DeepSeek/L5")
    p_dsr.add_argument("key", choices=["mdc-u", "libros-tabla", "criba-aleatorovix"])
    p_dsr.add_argument("--script", "-s", help="ID script (ej. mdc-v23, ksweep-py)")
    p_dsr.add_argument("--n", type=int, help="N para mdc-v23 vía -c")
    p_dsr.add_argument("extra", nargs="*", help="Args extra al script Python")
    p_dsr.set_defaults(func=cmd_deepseek_run)

    p_teo = sub.add_parser("teorema", help="Fichas teoremas/ (INDICE_MAESTRO)")
    p_teo_sub = p_teo.add_subparsers(dest="teorema_cmd", required=True)
    p_teo_sub.add_parser("list", help="Índice 01–30").set_defaults(func=cmd_teorema_list)
    p_ti = p_teo_sub.add_parser("info", help="Leer ficha teorema")
    p_ti.add_argument("numero", type=int, nargs="?", help="Número 1–30")
    p_ti.add_argument("--id", help="ID alternativo (ej. dos_primos)")
    p_ti.add_argument("--no-body", action="store_true", help="Solo metadatos")
    p_ti.set_defaults(func=cmd_teorema_info)

    p_game = sub.add_parser("game", help="MMO masivo — UDP estado + Grafcet + mundo")
    p_game_sub = p_game.add_subparsers(dest="game_cmd", required=True)
    p_gd = p_game_sub.add_parser(
        "demo",
        help="Servidor autoritativo + bots (B slot-ring → D Grafcet → WorldState)",
    )
    p_gd.add_argument("--port", type=int, default=3344)
    p_gd.add_argument("--duration", type=float, default=5.0)
    p_gd.add_argument("--players", type=int, default=128)
    p_gd.add_argument("--tick-hz", type=float, default=20.0)
    p_gd.add_argument("--batch-size", type=int, default=32)
    p_gd.add_argument("--shards", type=int, default=4, help="Particiones de mundo")
    p_gd.add_argument("--duplicate-ratio", type=float, default=0.15)
    p_gd.add_argument("--out", metavar="FILE", help="Exportar metricas JSON")
    p_gd.add_argument(
        "--no-bots",
        action="store_true",
        help="Solo servidor (esperar clientes UDP externos)",
    )
    p_gd.set_defaults(func=cmd_game_demo)

    p_gc = p_game_sub.add_parser(
        "cluster",
        help="MMO + hubs L3 WORK/RESULT por shard (cluster distribuido)",
    )
    p_gc.add_argument("--port", type=int, default=3344)
    p_gc.add_argument("--duration", type=float, default=5.0)
    p_gc.add_argument("--players", type=int, default=128)
    p_gc.add_argument("--tick-hz", type=float, default=20.0)
    p_gc.add_argument("--batch-size", type=int, default=32)
    p_gc.add_argument("--shards", type=int, default=4)
    p_gc.add_argument("--hubs", type=int, default=4, help="Hubs L3 (>= shards)")
    p_gc.add_argument("--duplicate-ratio", type=float, default=0.15)
    p_gc.add_argument("--out", metavar="FILE")
    p_gc.add_argument("--no-bots", action="store_true")
    p_gc.add_argument(
        "--hub-hosts",
        help="Hubs remotos CSV (ej. 192.168.1.10,192.168.1.11) — sin spawn local",
    )
    p_gc.add_argument("--hub-base-port", type=int, default=19701)
    p_gc.set_defaults(func=cmd_game_cluster)

    p_gb = p_game_sub.add_parser("bots", help="Lanzar bots contra servidor existente")
    p_gb.add_argument("--master", default="127.0.0.1")
    p_gb.add_argument("--port", type=int, default=3344)
    p_gb.add_argument("--players", type=int, default=64)
    p_gb.add_argument("--duration", type=float, default=5.0)
    p_gb.add_argument("--tick-hz", type=float, default=20.0)
    p_gb.add_argument("--duplicate-ratio", type=float, default=0.12)
    p_gb.set_defaults(func=cmd_game_bots)

    p_mech = sub.add_parser("mechanical", help="Regla mecanica toy")
    p_mech.add_argument("a", type=float)
    p_mech.add_argument("b", type=float)
    p_mech.add_argument("deltas", nargs="*", type=float, default=[0.1, 0.5, 1.0])
    p_mech.set_defaults(func=cmd_mechanical)

    p_native = sub.add_parser("native", help="Nucleo C unificado (antipc_native.dll)")
    p_native_sub = p_native.add_subparsers(dest="native_cmd", required=True)
    p_native_sub.add_parser("status", help="Backend y DLL cargada").set_defaults(
        func=cmd_native_status
    )
    p_nb = p_native_sub.add_parser("bench", help="Benchmark C vs Python")
    p_nb.add_argument("--limit", type=int, default=50_000)
    p_nb.add_argument("--mdc-n", type=int, default=1_047_029)
    p_nb.add_argument("--json", metavar="FILE")
    p_nb.set_defaults(func=cmd_native_bench)

    p_criba = sub.add_parser("criba", help="Criba Eratostenes (C nativo)")
    p_criba.add_argument("--limit", type=int, default=10_000)
    p_criba.add_argument("--hibrida", action="store_true", help="Criba hibrida VMA (C)")
    p_criba.add_argument("--modular6k", action="store_true", help="Criba modular 6k±1 (C)")
    p_criba.add_argument(
        "--desmemoriada",
        action="store_true",
        help="Criba desmemoriada VMA (vma-methods Python)",
    )
    p_criba.add_argument("--list", action="store_true", help="Mostrar lista de primos")
    p_criba.add_argument("--max-list", type=int, default=32)
    p_criba.set_defaults(func=cmd_criba)

    p_newton = sub.add_parser("newton", help="Newton Rapido + oraculo MEcuation (C)")
    p_newton.add_argument("E", type=float, help="Valor E (ej. 121)")
    p_newton.add_argument("--base", "-b", type=float, default=10.0)
    p_newton.add_argument(
        "--familia",
        "-f",
        default="cuadrados",
        choices=["general", "cuadrados", "cubos", "potencia", "kp", "mersenne"],
    )
    p_newton.add_argument("--n-exp", type=int, default=2, help="Exponente (familia potencia)")
    p_newton.add_argument("--k", type=float, default=1.0, help="Factor k (familia kp)")
    p_newton.set_defaults(func=cmd_newton)

    p_gm = sub.add_parser("geo-masivo", help="Aleatorovix + fase pi/e (C)")
    p_gm.add_argument("--text", "-t", help="Texto a cifrar/descifrar demo")
    p_gm.add_argument("--file", "-f", help="Archivo demo")
    p_gm.add_argument("--semilla", type=int, default=43210)
    p_gm.set_defaults(func=cmd_geo_masivo)

    p_geo = sub.add_parser("geo", help="Convergencia geometrica binaria (C)")
    p_geo.add_argument(
        "--demo",
        action="store_true",
        help="Demo Libro 4 (_demo_convergencia.py defaults)",
    )
    p_geo.add_argument("--bits", default="0101101011000010")
    p_geo.add_argument("--tales", default="3,5,8,13,21")
    p_geo.add_argument("--puntos", default="6,12,18")
    p_geo.set_defaults(func=cmd_geo)

    p_k3 = sub.add_parser("k3", help="HASHTOOLCODE K3 suite (hash, dedup, Grafcet)")
    p_k3_sub = p_k3.add_subparsers(dest="k3_cmd", required=True)

    p_k3_h = p_k3_sub.add_parser("hash", help="Hash buffer (--text o --file)")
    p_k3_h.add_argument("--text", "-t")
    p_k3_h.add_argument("--file", "-f")
    p_k3_h.add_argument("--seed", type=lambda x: int(x, 0), default=None)
    p_k3_h.set_defaults(func=cmd_k3_hash)

    p_k3_f = p_k3_sub.add_parser("file", help="Hash fichero streaming C")
    p_k3_f.add_argument("file")
    p_k3_f.add_argument("--seed", type=lambda x: int(x, 0), default=None)
    p_k3_f.add_argument(
        "--fingerprint",
        action="store_true",
        help="Tamaño + hash en una llamada (k3dedup)",
    )
    p_k3_f.set_defaults(func=cmd_k3_file)

    p_k3_sub.add_parser("verify", help="Paridad engine vs Python").set_defaults(
        func=cmd_k3_verify
    )

    p_k3_d = p_k3_sub.add_parser("dedup", help="Duplicados por tam+hash (k3dedup)")
    p_k3_d.add_argument("files", nargs="*", help="Rutas de archivos")
    p_k3_d.add_argument("--dir", "-d", help="Escanear directorio recursivo")
    p_k3_d.add_argument("--json", metavar="FILE")
    p_k3_d.set_defaults(func=cmd_k3_dedup)

    p_k3_r = p_k3_sub.add_parser("redundant", help="Hashes redundantes Grafcet")
    p_k3_r.add_argument("--text", "-t")
    p_k3_r.add_argument("--file", "-f")
    p_k3_r.add_argument("--replicas", type=int, default=3)
    p_k3_r.set_defaults(func=cmd_k3_redundant)

    p_k3_he = p_k3_sub.add_parser("heavy", help="Heavy hash 4 rondas Grafcet")
    p_k3_he.add_argument("--text", "-t")
    p_k3_he.add_argument("--file", "-f")
    p_k3_he.set_defaults(func=cmd_k3_heavy)

    p_k3_hm = p_k3_sub.add_parser("hamming", help="Distancia Hamming entre hashes hex")
    p_k3_hm.add_argument("hash_a")
    p_k3_hm.add_argument("hash_b")
    p_k3_hm.set_defaults(func=cmd_k3_hamming)

    p_k3_si = p_k3_sub.add_parser("simil", help="Similitud Jaccard por shingles K3")
    p_k3_si.add_argument("files", nargs="*", help="Archivos texto")
    p_k3_si.add_argument("--dir", "-d", help="Directorio .txt")
    p_k3_si.add_argument("--threshold", type=float, default=0.30)
    p_k3_si.set_defaults(func=cmd_k3_simil)

    p_k3_sr = p_k3_sub.add_parser("search", help="Índice invertido k3search")
    p_k3_sr_sub = p_k3_sr.add_subparsers(dest="search_cmd", required=True)
    p_k3_sri = p_k3_sr_sub.add_parser("index", help="Crear índice .k3idx")
    p_k3_sri.add_argument("--output", "-o", default="output/indice.k3idx")
    p_k3_sri.add_argument("files", nargs="*", help="Archivos a indexar")
    p_k3_sri.add_argument("--dir", "-d", help="Directorio recursivo")
    p_k3_sri.set_defaults(func=cmd_k3_search_index)
    p_k3_srq = p_k3_sr_sub.add_parser("query", help="Consultar índice")
    p_k3_srq.add_argument("index")
    p_k3_srq.add_argument("words", nargs="+")
    p_k3_srq.add_argument("--top", type=int, default=20)
    p_k3_srq.set_defaults(func=cmd_k3_search_query)

    p_k3_x = p_k3_sub.add_parser("stream-xor", help="Motor acordeon K3 XOR (33x1)")
    p_k3_x.add_argument("--text", "-t")
    p_k3_x.add_argument("--file", "-f")
    p_k3_x.add_argument("--base", type=int, default=33)
    p_k3_x.add_argument("--rel", type=int, default=1)
    p_k3_x.set_defaults(func=cmd_k3_stream_xor)

    p_ind = sub.add_parser("industrial", help="Inventario y auditoria industrial")
    p_ind_sub = p_ind.add_subparsers(dest="ind_cmd", required=True)

    p_inv = p_ind_sub.add_parser("inventory", help="Listar que hace y que tiene AntiPC")
    p_inv.add_argument("--out", metavar="FILE", help="Exportar inventario TXT")
    p_inv.add_argument("--json", metavar="FILE", help="Exportar inventario JSON")
    p_inv.add_argument("-q", "--quiet", action="store_true", help="Solo rutas export")
    p_inv.set_defaults(func=cmd_industrial_inventory)

    p_aud = p_ind_sub.add_parser("audit", help="Escanear directorio + telemetria + inventario")
    p_aud.add_argument("--dir", "-d", help="Directorio a auditar (default src/antipc)")
    p_aud.add_argument("--max-files", type=int, default=80)
    p_aud.add_argument("--workers", type=int, default=4)
    p_aud.add_argument("--out-dir", metavar="DIR", default="output")
    p_aud.add_argument("--full", action="store_true", help="Mostrar inventario completo")
    p_aud.set_defaults(func=cmd_industrial_audit)

    args = parser.parse_args(argv)
    if args.version and args.command is None:
        return cmd_version(args)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())