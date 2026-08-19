#!/usr/bin/env python3
"""Extract text from Apiñon.docx and all embedded docx files."""
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
W_T = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t"
W_TAB = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}tab"
W_BR = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}br"
W_P = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}p"


def extract_text_from_docx_bytes(data: bytes) -> str:
    try:
        import io
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            xml = zf.read("word/document.xml")
    except Exception as e:
        return f"[ERROR reading docx: {e}]"

    try:
        root = ET.fromstring(xml)
    except ET.ParseError as e:
        return f"[ERROR parsing XML: {e}]"

    paragraphs = []
    for p in root.iter(W_P):
        parts = []
        for el in p.iter():
            if el.tag == W_T and el.text:
                parts.append(el.text)
            elif el.tag == W_T and el.tail:
                parts.append(el.tail)
            elif el.tag == W_TAB:
                parts.append("\t")
            elif el.tag == W_BR:
                parts.append("\n")
        line = "".join(parts).strip()
        if line:
            paragraphs.append(line)
    return "\n".join(paragraphs)


def extract_text_from_docx_path(path: Path) -> str:
    with zipfile.ZipFile(path) as zf:
        xml = zf.read("word/document.xml")
    root = ET.fromstring(xml)
    paragraphs = []
    for p in root.iter(W_P):
        parts = []
        for el in p.iter():
            if el.tag == W_T and el.text:
                parts.append(el.text)
            elif el.tag == W_T and el.tail:
                parts.append(el.tail)
            elif el.tag == W_TAB:
                parts.append("\t")
            elif el.tag == W_BR:
                parts.append("\n")
        line = "".join(parts).strip()
        if line:
            paragraphs.append(line)
    return "\n".join(paragraphs)


def first_lines(text: str, n: int = 5, max_len: int = 120) -> str:
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()][:n]
    out = []
    for ln in lines:
        if len(ln) > max_len:
            ln = ln[: max_len - 3] + "..."
        out.append(ln)
    return " | ".join(out) if out else "(vacío)"


def guess_title(text: str) -> str:
    for ln in text.splitlines():
        ln = ln.strip()
        if len(ln) >= 8:
            return ln[:100]
    return "(sin título)"


def extract_ole_strings(path: Path, min_len: int = 20) -> list[str]:
    data = path.read_bytes()
    # UTF-16LE strings
    strings = set()
    i = 0
    while i < len(data) - min_len * 2:
        if data[i] >= 0x20 and data[i] < 0x7F and data[i + 1] == 0:
            start = i
            chars = []
            j = i
            while j < len(data) - 1:
                lo, hi = data[j], data[j + 1]
                if hi == 0 and 0x20 <= lo < 0x7F:
                    chars.append(chr(lo))
                    j += 2
                else:
                    break
            if len(chars) >= min_len:
                s = "".join(chars)
                if not s.startswith("Microsoft") and "OLE" not in s[:10]:
                    strings.add(s[:200])
            i = j
        else:
            i += 1

    # ASCII runs
    for m in re.finditer(rb"[\x20-\x7e]{30,}", data):
        s = m.group().decode("ascii", errors="ignore")
        if "Microsoft" not in s and "OLE" not in s:
            strings.add(s[:200])
    return sorted(strings)[:8]


def main():
    apinnon = Path(r"C:\Users\cuent\Desktop\archivos\Apiñon.docx")
    extract_dir = Path(r"C:\Users\cuent\AppData\Local\Temp\apinnon_extract")
    out_dir = Path(r"C:\Users\cuent\Desktop\teoremas\apinnon_extraido")
    out_dir.mkdir(parents=True, exist_ok=True)

    if not extract_dir.exists():
        import shutil
        extract_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(apinnon) as zf:
            zf.extractall(extract_dir)

    # Main document (from unpacked XML)
    main_xml = extract_dir / "word" / "document.xml"
    root = ET.fromstring(main_xml.read_bytes())
    paragraphs = []
    for p in root.iter(W_P):
        parts = []
        for el in p.iter():
            if el.tag == W_T and el.text:
                parts.append(el.text)
            elif el.tag == W_T and el.tail:
                parts.append(el.tail)
            elif el.tag == W_TAB:
                parts.append("\t")
            elif el.tag == W_BR:
                parts.append("\n")
        line = "".join(parts).strip()
        if line:
            paragraphs.append(line)
    main_text = "\n".join(paragraphs)

    (out_dir / "00_texto_principal.txt").write_text(main_text, encoding="utf-8")

    # Embeddings
    emb_dir = extract_dir / "word" / "embeddings"
    embeddings = sorted(emb_dir.glob("*.docx"), key=lambda p: p.name)
    index_lines = []
    index_lines.append("INDICE APINON.DOCX — extracción automática")
    index_lines.append(f"Fuente: {apinnon}")
    index_lines.append(f"Fecha extracción: 2026-07-10")
    index_lines.append(f"Texto principal: {len(main_text)} caracteres, {len(paragraphs)} párrafos")
    index_lines.append(f"Embeddings docx: {len(embeddings)}")
    index_lines.append("")

    for i, emb in enumerate(embeddings):
        text = extract_text_from_docx_path(emb)
        title = guess_title(text)
        fname = f"{i+1:02d}_{emb.stem}.txt"
        (out_dir / fname).write_text(text, encoding="utf-8")
        chars = len(text)
        paras = len([l for l in text.splitlines() if l.strip()])
        index_lines.append(f"--- {emb.name} ---")
        index_lines.append(f"  Archivo extraído: apinnon_extraido/{fname}")
        index_lines.append(f"  Título/heurística: {title}")
        index_lines.append(f"  Tamaño: {emb.stat().st_size} bytes | {chars} chars | {paras} párrafos")
        index_lines.append(f"  Inicio: {first_lines(text, 3)}")
        index_lines.append("")

    # OLE objects
    ole_dir = extract_dir / "word" / "embeddings"
    ole_bins = sorted((extract_dir / "word").glob("**/oleObject*.bin"), key=lambda p: p.name)
    if not ole_bins:
        ole_bins = sorted(extract_dir.rglob("oleObject*.bin"), key=lambda p: p.name)

    index_lines.append(f"=== OLE OBJECTS ({len(ole_bins)}) ===")
    ole_out = out_dir / "ole_objects"
    ole_out.mkdir(exist_ok=True)
    for ob in ole_bins:
        strings = extract_ole_strings(ob)
        preview = strings[0][:150] if strings else "(sin texto legible)"
        index_lines.append(f"--- {ob.name} ({ob.stat().st_size} bytes) ---")
        index_lines.append(f"  Preview: {preview}")
        if strings:
            (ole_out / f"{ob.stem}.txt").write_text("\n---\n".join(strings), encoding="utf-8")
        index_lines.append("")

    # Media count
    media = list((extract_dir / "word" / "media").glob("*")) if (extract_dir / "word" / "media").exists() else []
    index_lines.append(f"=== IMÁGENES word/media: {len(media)} archivos ===")

    index_path = Path(r"C:\Users\cuent\Desktop\teoremas\INDICE_APINON.txt")
    index_path.write_text("\n".join(index_lines), encoding="utf-8")
    print(f"Main: {len(main_text)} chars")
    print(f"Embeddings: {len(embeddings)}")
    print(f"OLE: {len(ole_bins)}")
    print(f"Media: {len(media)}")
    print(f"Index: {index_path}")


if __name__ == "__main__":
    main()