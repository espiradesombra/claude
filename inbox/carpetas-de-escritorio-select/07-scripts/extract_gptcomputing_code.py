#!/usr/bin/env python3
"""
Extrae bloques de codigo C/C++ del gptcomputing.txt a from_gptcomputing/cpp/.
Heuristica: class/struct/enum/#include/int main/while(runtime
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TXT = ROOT.parent / "gptcomputing" / "gptcomputing.txt"
OUT_CPP = ROOT / "from_gptcomputing" / "cpp"
OUT_IDX = ROOT / "from_gptcomputing" / "INDICE_CODIGO_TXT.txt"

START_RE = re.compile(
    r"^(#pragma once|#include |template<|class \w|struct \w|enum class \w|enum \w|int main\(|while\s*\()"
)
CODE_CONT = re.compile(
    r"^(\s*\{|\s*\}|\s*public:|\s*private:|\s*protected:|\s*virtual |\s*override|"
    r"\s*return |\s*if\s*\(|\s*for\s*\(|\s*using |\s*namespace |\s*//|"
    r"\s*[A-Za-z_][\w:<>,\s*&]*[;\{]|^\s*\};?\s*$)"
)

# Rangos manuales para bloques multi-seccion o headers completos (lineas 1-based del txt)
MANUAL_RANGES: dict[str, tuple[int, int]] = {
    "core/include/types/Reference.h": (1672, 1819),
    "core/src/storage/PayloadStore.cpp": (1911, 1950),
    "core/include/knowledge/KnowledgeRouter.h": (3381, 3415),
    "core/src/knowledge/KnowledgeRouter_resolve.cpp": (3430, 3475),
    "main.cpp": (2520, 2533),
    "k3/include/KOPHeader.h": (9802, 9809),
    "k3/include/KOP.h": (9754, 9758),
    "core/include/identity/KnowledgeDNA.h": (8436, 8495),
    "core/include/identity/Signature_class.h": (8294, 8365),
    "core/include/identity/KnowledgeID.h": (8195, 8231),
}


def is_code_line(line: str) -> bool:
    s = line.rstrip()
    if not s.strip():
        return True
    if START_RE.match(s.strip()):
        return True
    t = s.strip()
    if t in ("{", "}", "};", "public:", "private:", "protected:"):
        return True
    if any(
        x in t
        for x in (
            ";",
            "virtual ",
            "override",
            "namespace ",
            "#include",
            "return ",
            "std::",
            "uint",
            "bool ",
            "double ",
            "float ",
            "void ",
            "int ",
            "template",
            "constexpr",
            "operator",
            "noexcept",
            "/*",
            "//",
        )
    ):
        return True
    if re.match(r"^\s*[A-Za-z_][\w]*\s*[\(\{]", t):
        return True
    return False


def extract_auto_blocks(lines: list[str]) -> list[tuple[str, int, int, list[str]]]:
    blocks: list[tuple[str, int, int, list[str]]] = []
    i = 0
    n = len(lines)
    used_names: dict[str, int] = {}

    while i < n:
        raw = lines[i]
        stripped = raw.strip()
        if not START_RE.match(stripped):
            i += 1
            continue
        start = i + 1
        block = [raw]
        i += 1
        brace_depth = raw.count("{") - raw.count("}")
        # int main / while without class — short block
        if stripped.startswith("while") or stripped.startswith("int main"):
            while i < n and (is_code_line(lines[i]) or not lines[i].strip()):
                if lines[i].strip():
                    block.append(lines[i])
                    brace_depth += lines[i].count("{") - lines[i].count("}")
                i += 1
                if brace_depth <= 0 and len(block) > 2:
                    break
        else:
            while i < n:
                if not lines[i].strip():
                    if brace_depth <= 0 and len(block) >= 3:
                        break
                    block.append(lines[i])
                    i += 1
                    continue
                if not is_code_line(lines[i]) and brace_depth <= 0:
                    break
                block.append(lines[i])
                brace_depth += lines[i].count("{") - lines[i].count("}")
                i += 1
                if stripped.startswith("#include") and ";" in lines[i - 1]:
                    break

        content = "".join(block).strip()
        if len(content) < 20:
            continue
        name = _guess_filename(stripped, content)
        used_names[name] = used_names.get(name, 0) + 1
        if used_names[name] > 1:
            name = name.replace(".h", f"_{used_names[name]}.h")
        blocks.append((name, start, i, block))
    return blocks


def _guess_filename(first_line: str, content: str) -> str:
    m = re.search(r"class\s+(\w+)", content)
    if m:
        return f"core/include/auto/{m.group(1)}.h"
    m = re.search(r"struct\s+(\w+)", content)
    if m:
        return f"core/include/auto/{m.group(1)}.h"
    m = re.search(r"enum class\s+(\w+)", content)
    if m:
        return f"core/include/auto/{m.group(1)}.h"
    if first_line.startswith("int main"):
        return "main.cpp"
    if first_line.startswith("while"):
        return "core/include/runtime/Runtime_loop.h"
    if first_line.startswith("#include"):
        return "core/include/auto/includes_fragment.h"
    return "core/include/auto/fragment.h"


def wrap_header(path: str, body: str) -> str:
    if body.startswith("#pragma"):
        return body + "\n"
    guard = "ANTIPC_" + re.sub(r"[^A-Z0-9]", "_", path.upper()) + "_H"
    return (
        f"#pragma once\n"
        f"// Extraido de gptcomputing.txt -> {path}\n"
        f"// Implementacion Python activa: src/antipc/runtime/\n\n"
        f"#ifndef {guard}\n#define {guard}\n\n"
        f"namespace antipc {{\n\n"
        f"{body}\n\n"
        f"}} // namespace antipc\n\n"
        f"#endif\n"
    )


def main() -> None:
    lines = TXT.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    OUT_CPP.mkdir(parents=True, exist_ok=True)

    written: list[str] = []

    # Manual ranges first (canonical large blocks)
    for rel, (a, b) in MANUAL_RANGES.items():
        chunk = "".join(lines[a - 1 : b]).strip()
        path = OUT_CPP / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        if not chunk.startswith("#pragma"):
            chunk = wrap_header(rel, chunk)
        path.write_text(chunk + "\n", encoding="utf-8")
        written.append(f"L{a}-L{b}  {rel}")

    # Auto-extracted blocks
    auto = extract_auto_blocks(lines)
    for name, start, end, block in auto:
        rel = name
        full = OUT_CPP / rel
        if full.exists():
            rel = rel.replace(".h", f"_L{start}.h")
            full = OUT_CPP / rel
        full.parent.mkdir(parents=True, exist_ok=True)
        body = "".join(block).strip()
        if body.startswith("#include") and "class " not in body and "struct " not in body:
            text = (
                f"// Fragmento includes lineas {start}-{end}\n{body}\n"
            )
        else:
            text = wrap_header(rel, body)
        full.write_text(text, encoding="utf-8")
        written.append(f"L{start}-L{end}  {rel}")

    # Python map — codigo del txt implementado en runtime
    py_map = ROOT / "from_gptcomputing" / "PYTHON_IMPLEMENTADO.txt"
    py_map.write_text(
        """# gptcomputing.txt -> Python implementado (v0.1.0 / v0.1.1)

v0.0.01 Flow Kernel          src/antipc/runtime/kernel.py
v0.0.02 Knowledge Buffer     src/antipc/runtime/knowledge.py
v0.0.03 Reference Graph      src/antipc/runtime/reference.py
v0.0.04 Signature            src/antipc/runtime/signature.py
v0.0.065 IPlugin/Capability  src/antipc/runtime/plugin.py
v0.0.065 PluginManager       src/antipc/runtime/plugin_manager.py
v0.0.06 Modes                src/antipc/runtime/modes.py
EventBus                     src/antipc/runtime/event_bus.py
Scheduler                    src/antipc/runtime/scheduler.py
PendingOperation             src/antipc/runtime/scheduler.py
Resolver                     src/antipc/runtime/resolver.py
KOP 001-008 binario          src/antipc/runtime/kop.py
KOP Registry                 src/antipc/runtime/kop_registry.py
Arquitecturas A-E            src/antipc/architectures.py
UDP protocol                 src/antipc/udp_protocol.py
udp_benchmark                src/antipc/udp_benchmark.py
hub_node                     src/antipc/hub_node.py
HashPlugin/K3 plugins        src/antipc/plugins/
demo_v010                    src/antipc/demo_v010.py
demo_kop_binary              src/antipc/demo_kop_binary.py
""",
        encoding="utf-8",
    )
    written.append("PYTHON_IMPLEMENTADO.txt")

    idx_lines = [
        "=" * 72,
        "INDICE — codigo extraido de gptcomputing.txt",
        f"Origen: {TXT}",
        f"Destino: {OUT_CPP}",
        f"Total archivos: {len(written)}",
        "=" * 72,
        "",
        "C++ / pseudocodigo del chat (referencia Sprint 2+):",
    ]
    idx_lines.extend(f"  {w}" for w in sorted(written))
    idx_lines.extend(
        [
            "",
            "Python EJECUTABLE (no duplicado en cpp/):",
            "  Ver PYTHON_IMPLEMENTADO.txt",
            "",
            "Regenerar:",
            "  python scripts/extract_gptcomputing_code.py",
            "",
        ]
    )
    OUT_IDX.write_text("\n".join(idx_lines), encoding="utf-8")
    print(f"OK: {len(written)} entradas -> {OUT_CPP}")
    print(f"Indice: {OUT_IDX}")


if __name__ == "__main__":
    main()