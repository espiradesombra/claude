"""
HASHTOOLCODE k3simil + k3search — port Python (vma-k3/c/k3hash/examples/).
"""

from __future__ import annotations

import re
from pathlib import Path

from hash_engine import k3_hash_buffer

K_SHINGLE = 5
INDEX_MAGIC = "K3SEARCH_INDEX_V1"


def tokenize_words(text: str) -> list[str]:
    return re.findall(r"[a-z0-9]+", text.lower())


def read_text_file(path: str | Path) -> str:
    return Path(path).read_text(encoding="utf-8", errors="replace")


def build_shingle_hashes(text: str, k: int = K_SHINGLE) -> list[int]:
    words = tokenize_words(text)
    if len(words) < k:
        return []
    raw: list[int] = []
    for i in range(len(words) - k + 1):
        shingle = " ".join(words[i : i + k])
        raw.append(k3_hash_buffer(shingle.encode("utf-8")))
    return sorted(set(raw))


def jaccard(a: list[int], b: list[int]) -> float:
    if not a or not b:
        return 0.0
    i = j = inter = 0
    while i < len(a) and j < len(b):
        if a[i] == b[j]:
            inter += 1
            i += 1
            j += 1
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1
    union = len(a) + len(b) - inter
    return inter / union if union else 0.0


def k3_simil_pairs(
    paths: list[str | Path],
    *,
    threshold: float = 0.30,
    k: int = K_SHINGLE,
) -> list[dict]:
    docs: list[tuple[str, list[int]]] = []
    for raw in paths:
        p = Path(raw)
        if not p.is_file():
            continue
        shingles = build_shingle_hashes(read_text_file(p), k=k)
        docs.append((str(p.resolve()), shingles))

    pairs: list[dict] = []
    for i in range(len(docs)):
        for j in range(i + 1, len(docs)):
            sim = jaccard(docs[i][1], docs[j][1])
            if sim >= threshold:
                pairs.append(
                    {
                        "similarity": round(sim, 4),
                        "path_a": docs[i][0],
                        "path_b": docs[j][0],
                    }
                )
    pairs.sort(key=lambda x: -x["similarity"])
    return pairs


def k3_search_index(paths: list[str | Path], out_path: str | Path) -> dict:
    """Índice invertido texto plano K3SEARCH_INDEX_V1."""
    rutas: list[str] = []
    term_map: dict[str, dict[int, int]] = {}

    for raw in paths:
        p = Path(raw)
        if not p.is_file():
            continue
        doc_id = len(rutas)
        rutas.append(str(p.resolve()))
        for word in tokenize_words(read_text_file(p)):
            term_map.setdefault(word, {})
            term_map[word][doc_id] = term_map[word].get(doc_id, 0) + 1

    lines = [INDEX_MAGIC, f"DOCS {len(rutas)}"]
    for i, ruta in enumerate(rutas):
        lines.append(f"{i}\t{ruta}")
    lines.append(f"TERMS {len(term_map)}")
    for term in sorted(term_map.keys()):
        postings = term_map[term]
        parts = [term, str(len(postings))]
        for doc_id, cnt in sorted(postings.items()):
            parts.append(f"{doc_id}:{cnt}")
        lines.append("\t".join(parts))

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"documents": len(rutas), "terms": len(term_map), "index": str(out)}


def _load_index(path: Path) -> tuple[list[str], dict[str, dict[int, int]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != INDEX_MAGIC:
        raise ValueError(f"índice inválido: {path}")
    idx = 1
    num_docs = int(lines[idx].split()[1])
    idx += 1
    rutas: list[str] = [""] * num_docs
    for _ in range(num_docs):
        doc_id_s, ruta = lines[idx].split("\t", 1)
        rutas[int(doc_id_s)] = ruta
        idx += 1
    num_terms = int(lines[idx].split()[1])
    idx += 1
    term_map: dict[str, dict[int, int]] = {}
    for _ in range(num_terms):
        parts = lines[idx].split("\t")
        term = parts[0]
        npost = int(parts[1])
        postings: dict[int, int] = {}
        for p in parts[2 : 2 + npost]:
            doc_s, cnt_s = p.split(":")
            postings[int(doc_s)] = int(cnt_s)
        term_map[term] = postings
        idx += 1
    return rutas, term_map


def k3_search_query(index_path: str | Path, words: list[str], *, top: int = 20) -> list[dict]:
    rutas, term_map = _load_index(Path(index_path))
    scores = [0.0] * len(rutas)
    for raw in words:
        term = raw.lower()
        if term not in term_map:
            continue
        for doc_id, cnt in term_map[term].items():
            scores[doc_id] += cnt
    ranked = sorted(
        [(i, scores[i]) for i in range(len(rutas)) if scores[i] > 0],
        key=lambda x: -x[1],
    )
    return [
        {"score": sc, "path": rutas[i]}
        for i, sc in ranked[:top]
    ]