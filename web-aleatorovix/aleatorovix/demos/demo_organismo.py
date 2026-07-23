#!/usr/bin/env python3
"""Demo: 10 ciclos como organismo_lila_v99.c"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from aleatorovix import AleatorovixOrganismo

if __name__ == "__main__":
    AleatorovixOrganismo().respirar(n=10)