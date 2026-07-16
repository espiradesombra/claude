"""Telemetría industrial — auditoría CSV compatible SCADA/laboratorio."""

from __future__ import annotations

import csv
import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass
class TelemetryRecord:
    timestamp: float
    event: str
    plugin: str
    signature: str
    status: str
    elapsed_ms: float
    alu_saved: bool


@dataclass
class IndustrialTelemetry:
    records: list[TelemetryRecord] = field(default_factory=list)
    session_id: str = field(default_factory=lambda: str(int(time.time())))

    def log(self, event: str, plugin: str, signature: str, status: str,
            elapsed_ms: float = 0.0, alu_saved: bool = False) -> None:
        self.records.append(TelemetryRecord(
            time.time(), event, plugin, signature[:16], status, elapsed_ms, alu_saved
        ))

    def export_csv(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["timestamp", "event", "plugin", "signature", "status",
                        "elapsed_ms", "alu_saved"])
            for r in self.records:
                w.writerow([r.timestamp, r.event, r.plugin, r.signature,
                            r.status, f"{r.elapsed_ms:.3f}", int(r.alu_saved)])

    def export_json(self, path: Path, extra: dict | None = None) -> None:
        payload = {
            "session_id": self.session_id,
            "records": [asdict(r) for r in self.records],
            "summary": extra or {},
        }
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")