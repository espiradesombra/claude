#!/usr/bin/env python3
"""
EliminaDuplicadosK3 — GUI portable (hash K3).
Escanea una carpeta, agrupa duplicados exactos (tamaño + K3) y solo borra
tras confirmación explícita del usuario. Siempre deja al menos 1 copia.
"""

from __future__ import annotations

import os
import sys
import threading
import tkinter as tk
from collections import defaultdict
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

# --- K3 (mismo algoritmo que k3hash.c / antipc) ---
OFFSETS = (5, 7, 9, 11, 13, 17, 19, 23)
MAGICO = 0x9E3779B1
FINAL_MUL_A = 0x85EBCA6B
FINAL_MUL_B = 0xC2B2AE35
SEED = 0x1F2E3D4C
N_REGS = 4
ANCHO = 4

SKIP_DIR_NAMES = {".git", "node_modules", "__pycache__", ".venv", "venv", "$RECYCLE.BIN", "System Volume Information"}


def _rotl32(value: int, positions: int) -> int:
    positions &= 31
    value &= 0xFFFFFFFF
    if positions == 0:
        return value
    return ((value << positions) | (value >> (32 - positions))) & 0xFFFFFFFF


def _compress(estado: list[int], b: int) -> None:
    n = N_REGS
    b &= 0xFFFFFFFF
    x = (estado[0] ^ (b & estado[1 % n])) & 0xFFFFFFFF
    for j in range(2, n):
        x ^= _rotl32(estado[j], OFFSETS[j % 8])
    x ^= _rotl32(estado[1 % n], OFFSETS[4])
    x = (x + _rotl32(estado[0], OFFSETS[5])) & 0xFFFFFFFF
    x ^= (b * MAGICO) & 0xFFFFFFFF
    x ^= b
    x = (x + _rotl32(b, OFFSETS[6])) & 0xFFFFFFFF
    x ^= _rotl32(b, OFFSETS[7])
    x ^= _rotl32(x, OFFSETS[2])
    x = (x * MAGICO) & 0xFFFFFFFF
    x ^= x >> 15
    prev0 = estado[0]
    estado[0] = x & 0xFFFFFFFF
    temp_prev = prev0
    for i in range(1, n):
        temp_actual = estado[i]
        estado[i] = (_rotl32(estado[i], OFFSETS[0] + i) ^ temp_prev) & 0xFFFFFFFF
        temp_prev = temp_actual
        if i == n - 1:
            estado[i] = (estado[i] ^ b) & 0xFFFFFFFF


def _finalize(regs: list[int]) -> int:
    acc = regs[0]
    for i in range(1, N_REGS):
        acc ^= _rotl32(regs[i], 5 + i)
        acc &= 0xFFFFFFFF
    acc ^= acc >> 16
    acc = (acc * FINAL_MUL_A) & 0xFFFFFFFF
    acc ^= acc >> 13
    acc = (acc * FINAL_MUL_B) & 0xFFFFFFFF
    acc ^= acc >> 16
    return acc & 0xFFFFFFFF


def k3_hash_file(path: Path, chunk: int = 1 << 20) -> int:
    regs = [(SEED ^ (i * MAGICO)) & 0xFFFFFFFF for i in range(N_REGS)]
    partial = bytearray()
    with path.open("rb") as f:
        while True:
            data = f.read(chunk)
            if not data:
                break
            buf = partial + data
            i = 0
            while i + ANCHO <= len(buf):
                bloque = 0
                for b in range(ANCHO):
                    bloque = (bloque << 8) | buf[i + b]
                _compress(regs, bloque)
                i += ANCHO
            partial = bytearray(buf[i:])
    if partial:
        bloque = 0
        for b in range(ANCHO):
            byte = partial[b] if b < len(partial) else 0
            bloque = (bloque << 8) | byte
        _compress(regs, bloque)
    return _finalize(regs)


def human_bytes(n: int) -> str:
    x = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(x) < 1024.0:
            return f"{x:.1f} {unit}"
        x /= 1024.0
    return f"{x:.1f} PB"


def collect_files(root: Path, progress_cb=None) -> list[Path]:
    out: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIR_NAMES and not d.startswith(".")]
        base = Path(dirpath)
        for name in filenames:
            out.append(base / name)
            if progress_cb and len(out) % 200 == 0:
                progress_cb(f"Listando… {len(out)} ficheros")
    return out


def find_duplicate_groups(files: list[Path], progress_cb=None) -> list[dict]:
    by_size: dict[int, list[Path]] = defaultdict(list)
    for p in files:
        try:
            by_size[p.stat().st_size].append(p)
        except OSError:
            continue

    groups: list[dict] = []
    candidates = [(s, ps) for s, ps in by_size.items() if len(ps) > 1]
    hashed = 0
    total_hash = sum(len(ps) for _, ps in candidates)

    for size, paths in candidates:
        by_hash: dict[int, list[Path]] = defaultdict(list)
        for p in paths:
            try:
                h = k3_hash_file(p)
                hashed += 1
                by_hash[h].append(p)
                if progress_cb and hashed % 5 == 0:
                    progress_cb(f"Hasheando K3… {hashed}/{total_hash}")
            except OSError:
                continue
        for h, same in by_hash.items():
            if len(same) < 2:
                continue
            # Conservar por defecto: ruta más corta, luego más antigua mtime
            ranked = sorted(
                same,
                key=lambda x: (
                    len(str(x)),
                    x.stat().st_mtime if x.exists() else 0,
                    str(x).lower(),
                ),
            )
            keep = ranked[0]
            groups.append(
                {
                    "size": size,
                    "hash": h,
                    "paths": ranked,
                    "keep": keep,
                    "delete_defaults": ranked[1:],
                }
            )

    groups.sort(key=lambda g: -(g["size"] * (len(g["paths"]) - 1)))
    return groups


class App(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("EliminaDuplicados K3 — VMA / HASHTOOLCODE")
        self.geometry("980x640")
        self.minsize(720, 480)

        self.folder = tk.StringVar(value="")
        self.status = tk.StringVar(value="Elige una carpeta y pulsa Escanear.")
        self.groups: list[dict] = []
        # iid -> (group_idx, path or None if header)
        self._row_meta: dict[str, tuple[int, Path | None]] = {}
        self._scanning = False

        top = ttk.Frame(self, padding=8)
        top.pack(fill=tk.X)
        ttk.Label(top, text="Carpeta:").pack(side=tk.LEFT)
        ttk.Entry(top, textvariable=self.folder).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=6)
        ttk.Button(top, text="Examinar…", command=self.pick_folder).pack(side=tk.LEFT)
        ttk.Button(top, text="Escanear", command=self.start_scan).pack(side=tk.LEFT, padx=(6, 0))

        mid = ttk.Frame(self, padding=(8, 0))
        mid.pack(fill=tk.BOTH, expand=True)

        cols = ("mark", "role", "size", "hash", "path")
        self.tree = ttk.Treeview(mid, columns=cols, show="headings", selectmode="extended")
        self.tree.heading("mark", text="Borrar?")
        self.tree.heading("role", text="Rol")
        self.tree.heading("size", text="Tamaño")
        self.tree.heading("hash", text="Hash K3")
        self.tree.heading("path", text="Ruta")
        self.tree.column("mark", width=70, anchor=tk.CENTER)
        self.tree.column("role", width=90, anchor=tk.CENTER)
        self.tree.column("size", width=90, anchor=tk.E)
        self.tree.column("hash", width=110, anchor=tk.CENTER)
        self.tree.column("path", width=560, anchor=tk.W)
        vsb = ttk.Scrollbar(mid, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.bind("<Double-1>", self.toggle_mark)
        self.tree.bind("<space>", self.toggle_mark)

        hint = ttk.Label(
            self,
            text="Doble clic o Espacio: marcar/desmarcar borrado. En cada grupo siempre debe quedar ≥1 archivo (CONSERVAR).",
            padding=(8, 4),
        )
        hint.pack(fill=tk.X)

        bot = ttk.Frame(self, padding=8)
        bot.pack(fill=tk.X)
        ttk.Button(bot, text="Marcar todos los duplicados", command=self.mark_all_dupes).pack(side=tk.LEFT)
        ttk.Button(bot, text="Desmarcar todo", command=self.unmark_all).pack(side=tk.LEFT, padx=6)
        ttk.Button(bot, text="⚠ Eliminar marcados…", command=self.delete_marked).pack(side=tk.RIGHT)
        ttk.Label(self, textvariable=self.status, padding=(8, 4)).pack(fill=tk.X)

    def pick_folder(self) -> None:
        path = filedialog.askdirectory(title="Carpeta a escanear")
        if path:
            self.folder.set(path)

    def start_scan(self) -> None:
        if self._scanning:
            return
        root = Path(self.folder.get().strip())
        if not root.is_dir():
            messagebox.showerror("Error", "Elige una carpeta válida.")
            return
        self._scanning = True
        self.tree.delete(*self.tree.get_children())
        self.groups = []
        self._row_meta.clear()
        self.status.set("Escaneando…")

        def work() -> None:
            try:
                files = collect_files(root, lambda m: self.after(0, lambda: self.status.set(m)))
                groups = find_duplicate_groups(files, lambda m: self.after(0, lambda: self.status.set(m)))
                self.after(0, lambda: self._load_groups(groups, len(files)))
            except Exception as e:
                self.after(0, lambda: messagebox.showerror("Error", str(e)))
            finally:
                self.after(0, lambda: setattr(self, "_scanning", False))

        threading.Thread(target=work, daemon=True).start()

    def _load_groups(self, groups: list[dict], nfiles: int) -> None:
        self.groups = groups
        recoverable = sum(g["size"] * (len(g["paths"]) - 1) for g in groups)
        self.status.set(
            f"{nfiles} ficheros · {len(groups)} grupos de duplicados · ~{human_bytes(recoverable)} recuperables"
        )
        for gi, g in enumerate(groups):
            header = self.tree.insert(
                "",
                tk.END,
                values=(
                    "",
                    f"GRUPO {gi + 1}",
                    human_bytes(g["size"]),
                    f"0x{g['hash']:08X}",
                    f"{len(g['paths'])} copias — clic para expandir selección",
                ),
                tags=("header",),
            )
            self._row_meta[header] = (gi, None)
            for p in g["paths"]:
                is_keep = p == g["keep"]
                mark = "" if is_keep else "☑"
                role = "CONSERVAR" if is_keep else "borrar?"
                iid = self.tree.insert(
                    "",
                    tk.END,
                    values=(mark, role, human_bytes(g["size"]), f"0x{g['hash']:08X}", str(p)),
                    tags=("keep",) if is_keep else ("dupe",),
                )
                self._row_meta[iid] = (gi, p)
        self.tree.tag_configure("header", background="#e8e8e8")
        self.tree.tag_configure("keep", background="#e6ffe6")
        self.tree.tag_configure("dupe", background="#fff5e6")

    def toggle_mark(self, _evt=None) -> None:
        sel = self.tree.selection()
        if not sel:
            return
        for iid in sel:
            meta = self._row_meta.get(iid)
            if not meta or meta[1] is None:
                continue
            gi, path = meta
            g = self.groups[gi]
            vals = list(self.tree.item(iid, "values"))
            currently = vals[0] == "☑"
            if currently:
                vals[0] = ""
                vals[1] = "CONSERVAR" if path == g["keep"] else "mantener"
                self.tree.item(iid, values=vals, tags=("keep",))
            else:
                # No permitir marcar si sería el único no marcado del grupo
                others_marked_or_keep = 0
                for other_iid, (ogi, op) in self._row_meta.items():
                    if ogi != gi or op is None or other_iid == iid:
                        continue
                    ov = self.tree.item(other_iid, "values")
                    if ov[0] != "☑":
                        others_marked_or_keep += 1
                if others_marked_or_keep < 1:
                    messagebox.showwarning(
                        "Atención",
                        "No puedes marcar para borrar todas las copias de un grupo.\n"
                        "Debe quedar al menos una.",
                    )
                    return
                vals[0] = "☑"
                vals[1] = "borrar?"
                self.tree.item(iid, values=vals, tags=("dupe",))

    def mark_all_dupes(self) -> None:
        for iid, (gi, path) in list(self._row_meta.items()):
            if path is None:
                continue
            g = self.groups[gi]
            if path == g["keep"]:
                continue
            vals = list(self.tree.item(iid, "values"))
            vals[0] = "☑"
            vals[1] = "borrar?"
            self.tree.item(iid, values=vals, tags=("dupe",))

    def unmark_all(self) -> None:
        for iid, (gi, path) in list(self._row_meta.items()):
            if path is None:
                continue
            g = self.groups[gi]
            vals = list(self.tree.item(iid, "values"))
            vals[0] = ""
            vals[1] = "CONSERVAR" if path == g["keep"] else "mantener"
            self.tree.item(iid, values=vals, tags=("keep",))

    def _marked_paths(self) -> list[Path]:
        out: list[Path] = []
        for iid, (_gi, path) in self._row_meta.items():
            if path is None:
                continue
            vals = self.tree.item(iid, "values")
            if vals[0] == "☑":
                out.append(path)
        return out

    def delete_marked(self) -> None:
        targets = self._marked_paths()
        if not targets:
            messagebox.showinfo("Nada que borrar", "No hay ficheros marcados con ☑.")
            return

        # Safety: per group at least one survivor
        for gi, g in enumerate(self.groups):
            survivors = []
            doomed = []
            for iid, (ogi, path) in self._row_meta.items():
                if ogi != gi or path is None:
                    continue
                if self.tree.item(iid, "values")[0] == "☑":
                    doomed.append(path)
                else:
                    survivors.append(path)
            if doomed and not survivors:
                messagebox.showerror(
                    "Bloqueado",
                    f"El grupo {gi + 1} borraría todas las copias. Desmarca al menos una.",
                )
                return

        total = sum(p.stat().st_size for p in targets if p.exists())
        preview = "\n".join(str(p) for p in targets[:15])
        if len(targets) > 15:
            preview += f"\n… y {len(targets) - 15} más"

        if not messagebox.askokcancel(
            "Confirmación 1/2",
            f"Vas a borrar {len(targets)} fichero(s) (~{human_bytes(total)}).\n\n"
            f"{preview}\n\n"
            "Esto NO se puede deshacer desde esta app.\n¿Continuar a la confirmación final?",
        ):
            return

        # Confirmación fuerte: ventana que exige escribir BORRAR
        if not self._confirm_type_borrar(len(targets), total):
            messagebox.showinfo("Cancelado", "No se ha borrado nada.")
            return

        deleted = 0
        errors: list[str] = []
        for p in targets:
            try:
                p.unlink()
                deleted += 1
            except OSError as e:
                errors.append(f"{p}: {e}")

        msg = f"Borrados: {deleted} de {len(targets)}."
        if errors:
            msg += "\n\nErrores:\n" + "\n".join(errors[:10])
            messagebox.showwarning("Hecho con avisos", msg)
        else:
            messagebox.showinfo("Hecho", msg)

        # Rescan
        self.start_scan()

    def _confirm_type_borrar(self, n: int, total: int) -> bool:
        dlg = tk.Toplevel(self)
        dlg.title("Confirmación final — escribe BORRAR")
        dlg.transient(self)
        dlg.grab_set()
        dlg.geometry("480x220")
        result = {"ok": False}

        ttk.Label(
            dlg,
            text=(
                f"Para borrar {n} ficheros (~{human_bytes(total)}),\n"
                "escribe exactamente BORRAR (mayúsculas) y pulsa Confirmar.\n\n"
                "Cierra esta ventana o pulsa Cancelar para abortar."
            ),
            padding=12,
            justify=tk.CENTER,
        ).pack(fill=tk.X)

        entry = ttk.Entry(dlg, font=("Segoe UI", 14))
        entry.pack(padx=20, fill=tk.X)
        entry.focus_set()

        btns = ttk.Frame(dlg, padding=12)
        btns.pack(fill=tk.X)

        def do_ok() -> None:
            if entry.get().strip() == "BORRAR":
                result["ok"] = True
                dlg.destroy()
            else:
                messagebox.showerror("Incorrecto", 'Debes escribir exactamente "BORRAR".', parent=dlg)

        def do_cancel() -> None:
            dlg.destroy()

        ttk.Button(btns, text="Cancelar", command=do_cancel).pack(side=tk.LEFT)
        ttk.Button(btns, text="Confirmar borrado", command=do_ok).pack(side=tk.RIGHT)
        entry.bind("<Return>", lambda _e: do_ok())
        self.wait_window(dlg)
        return result["ok"]


def main() -> int:
    app = App()
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
