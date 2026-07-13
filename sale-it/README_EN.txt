================================================================================
  AntiPC — Adaptive Network Through Parallel Computing
  Sale-It Package v0.1
  Author: Víctor Manzanares Alberola
================================================================================

PACKAGE CONTENTS
--------------------------------------------------------------------------------

  src/antipc/          Python source code (full engine)
  docs/                Scientific documentation and benchmark results
  referencia/          Original K3 header (hash provenance)
  scripts/             Quick install and run scripts (Windows)
  README_EN.txt        This file (English)
  LEEME.txt            Spanish readme
  PRODUCT_SHEET_EN.txt Product sheet for international partners
  FICHA_PRODUCTO.txt   Product sheet (Spanish)
  LICENCIA.txt         License terms (Spanish)


REQUIREMENTS
--------------------------------------------------------------------------------

  - Windows 10/11 (or Linux/macOS with Python 3.10+)
  - Python 3.10 or newer
  - No external dependencies (standard library only)


QUICK START (WINDOWS)
--------------------------------------------------------------------------------

  1. Double-click:  INICIO.bat
  2. Choose option 1 to verify Python
  3. Choose option 3 for the recommended UDP demo (architectures A→E)


MANUAL COMMANDS
--------------------------------------------------------------------------------

  cd src\antipc
  python benchmark.py
  python benchmark.py --udp --packets 2000 --hubs 4
  python udp_benchmark.py --packets 2000 --hubs 4
  python verify_k3.py


ARCHITECTURES
--------------------------------------------------------------------------------

  A  Conventional     memcpy + blocking queue
  B  Lock-free        SPSC ring buffer, no extra copy
  C  Distributed      hub partitioning (threads) + K(N) cache
  D  Grafcet          batches + existence matrix + redundancy
  E  Real UDP         network hubs as compute nodes (localhost or LAN)


AntiPC FORMULA
--------------------------------------------------------------------------------

  P_util(N) = N · E(N) + K(N)

  N      = compute nodes / hubs
  E(N)   = coordination efficiency (0..1)
  K(N)   = reusable derived knowledge (cache hits)

================================================================================