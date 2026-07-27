#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
describe_folders.py — documenta el monorepo claude (espiradesombra/VMA).

Genera / actualiza:
  - MAPA.md          mapa temático (empezar aquí)
  - FOLDERS.md       índice alfabético + conteos
  - <carpeta>/README.md  descripción de carpeta y contenido (si es auto o --force)

No pisa READMEs "propios" (largos / manuales) salvo --force.
No lista chats personales ni teléfonos en los READMEs.

Uso:
  python tools/describe_folders.py
  python tools/describe_folders.py --force-weak
  python tools/describe_folders.py --only "junto,xd,Nueva carpeta"
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO_URL = "https://github.com/espiradesombra/claude"

# --- descripciones curadas (nombre exacto de carpeta de primer nivel) ---
CURATED: dict[str, dict] = {
    "1": {
        "title": "1 — núcleo personal / hilo 1477",
        "blurb": "Material canónico del autor: Quijote·Kilómetro, demos cuánticas toy, papers 3vs7, teorema PhaseAmplifier K3 XOR y chats de trabajo.",
        "tema": "Autor / corpus núcleo",
    },
    "2026": {
        "title": "2026 — línea de trabajo del año",
        "blurb": "Gemelos ZypyZape v3–v94, contexto Quijote, papers bilingües y capturas de control energético fechados en 2026.",
        "tema": "Gemelos / energía / 2026",
    },
    "33x1": {
        "title": "33×1 — el trato",
        "blurb": "Definición del trato: el **1** es TODO este monorepo civil; el **33** son 33 años de paz firmada por los países. Uso civil, comandos, manifiesto de hashes.",
        "tema": "33×1 (prioridad)",
    },
    "33x1-pack": {
        "title": "33×1-pack — paquete de distribución",
        "blurb": "Pack empaquetado del marco 33×1 para compartir / firmar (documentos de firma, material selectivo).",
        "tema": "33×1",
    },
    "33x1_CIENCIA": {
        "title": "33×1_CIENCIA — capa científica",
        "blurb": "Capa científica / anexos del objetivo 33×1.",
        "tema": "33×1",
    },
    "VMA_mates_rescat_2026": {
        "title": "VMA_mates_rescat_2026 — pack mates rescatado",
        "blurb": "Pack limpio 2026-07-23: sunraman, cribas (fix modular 6k±1), Criva/MRAUV, MDC, Newton Rápido, Sofí/Fermat/Goldbach, UPV C, gemelos y XFI N=3/4.",
        "tema": "Métodos / mates",
    },
    "aleatorovix": {
        "title": "aleatorovix — organismo de entropía",
        "blurb": "Organismo Aleatorovix: entropía, máscara Lila, criba desmemoriada, MDC, demos y benchmarks. Base del motor web sin Math.random.",
        "tema": "Aleatorovix",
    },
    "web-aleatorovix": {
        "title": "web-aleatorovix — ramo en el navegador",
        "blurb": "Sitio del ramo ❤️🌹 + motor JS Aleatorovix para techamv.com (sin Math.random). Incluye HTML, favicons y paquete técnico.",
        "tema": "Aleatorovix / web",
    },
    "techamv-aleatorovix-DEPLOY": {
        "title": "techamv-aleatorovix-DEPLOY — pack FTP mínimo",
        "blurb": "Copia mínima para subir por FTP a techamv.com: index.html + aleatorovix/aleatorovix.js (+ iconos).",
        "tema": "Web / deploy",
    },
    "techamv-web": {
        "title": "techamv-web — sitio ZypyZape / TechAMV",
        "blurb": "Sitio histórico ZypyZape/TechAMV: índices HTML, factor.php, chat y docs VMA.",
        "tema": "Web / ZypyZape",
    },
    "webtechamv": {
        "title": "webtechamv — hubs K3 / promo",
        "blurb": "Hubs K3, licencia software (React/TSX) y materiales promo HashCode.",
        "tema": "Web / K3",
    },
    "antipc": {
        "title": "antipc — Adaptive Network Through Parallel Computing",
        "blurb": "Motor industrial AntiPC: CLI, red UDP, MMO, DLL ALU, docs de arquitectura. Fórmula P_util(N)=N·E(N)+K(N).",
        "tema": "AntiPC",
    },
    "antipc-port-c": {
        "title": "antipc-port-c — port C de módulos",
        "blurb": "Port a C de criba, MDC trenes, Newton rápido + inventarios de integración v04–v14.",
        "tema": "AntiPC / C",
    },
    "antipc2": {
        "title": "antipc2 — línea secundaria AntiPC",
        "blurb": "Evolución / experimentos paralelos al árbol principal antipc/.",
        "tema": "AntiPC",
    },
    "cpu-antipc": {
        "title": "cpu-antipc — runtime CPU",
        "blurb": "Enfoque CPU/rendimiento y runtime AntiPC (concepto + C++).",
        "tema": "AntiPC",
    },
    "encriptacionGeometrica": {
        "title": "encriptacionGeometrica — cifrado por convergencia",
        "blurb": "Encriptación por convergencia geométrica (Thales, perímetros, K3 base 33 rel 1), demos Python/C y launcher.",
        "tema": "K3 / cifrado",
    },
    "vma-k3": {
        "title": "vma-k3 — motor K3",
        "blurb": "Integración VMA × K3: hash, dedup, demos C/Python, telemetría de marcas de código.",
        "tema": "K3 / cifrado",
    },
    "hashtool-work": {
        "title": "hashtool-work — taller de hash tools",
        "blurb": "Área de trabajo de herramientas hash / K3 (headers, C, notas).",
        "tema": "K3 / cifrado",
    },
    "hashtool-extract": {
        "title": "hashtool-extract — extractos hash/K3",
        "blurb": "Extracción y material de referencia de herramientas hash / K3.",
        "tema": "K3 / cifrado",
    },
    "vma-methods": {
        "title": "vma-methods — librería de métodos",
        "blurb": "Librería ejecutable: cribas (Desmemoriada, Modular 6k±1, Híbrida), Criva, Newton Rápido con oráculos.",
        "tema": "Métodos / mates",
    },
    "Metodo Newton Rápido": {
        "title": "Metodo Newton Rápido",
        "blurb": "Newton Rápido: logs, raíces, oráculos MEcuation, comparativas y código C++/Python.",
        "tema": "Métodos / mates",
    },
    "Metodo densidad MRAUV": {
        "title": "Metodo densidad MRAUV",
        "blurb": "Modelo de densidad de primos MRAUV (3 puntos cinemáticos) y documentos.",
        "tema": "Métodos / mates",
    },
    "mrauv": {
        "title": "mrauv — datos y scripts densitat",
        "blurb": "CSV y scripts del modelo MRAUV de densidad de primos.",
        "tema": "Métodos / mates",
    },
    "teoremas": {
        "title": "teoremas — colección formal VMA",
        "blurb": "Fichas de teoremas, cotas y resultados (formal / espacial / ingeniería).",
        "tema": "Métodos / mates",
    },
    "teoremasgrok": {
        "title": "teoremasgrok — teoremas con Grok",
        "blurb": "Derivaciones y fichas trabajadas con Grok (txt + headers C++).",
        "tema": "Métodos / mates",
    },
    "Goldbach": {
        "title": "Goldbach",
        "blurb": "Trabajo sobre la conjetura de Goldbach y variantes simétricas (docs, ppt).",
        "tema": "Métodos / mates",
    },
    "Sophie Germain": {
        "title": "Sophie Germain",
        "blurb": "Notas y presentaciones sobre Sophie Germain y estructura L₁ / 6k±1.",
        "tema": "Métodos / mates",
    },
    "sofi": {
        "title": "sofi — notas cortas Sophie",
        "blurb": "Material corto / LEEME de la línea Sophie / sofi structure.",
        "tema": "Métodos / mates",
    },
    "dos_primos": {
        "title": "dos_primos",
        "blurb": "Problemas de dos primos / primos gemelos (LEEME).",
        "tema": "Métodos / mates",
    },
    "siguiente_primo": {
        "title": "siguiente_primo",
        "blurb": "Búsqueda / validación del siguiente primo (fixes de paridad y fase).",
        "tema": "Métodos / mates",
    },
    "salto": {
        "title": "salto — saltos de primos",
        "blurb": "Saltos entre primos / 6k±1 (LEEME).",
        "tema": "Métodos / mates",
    },
    "discriminant": {
        "title": "discriminant",
        "blurb": "Demo de factorización por discriminantes y estructuras algebraicas.",
        "tema": "Métodos / mates",
    },
    "densidad": {
        "title": "densidad",
        "blurb": "Estudios y documentos de densidad de primos.",
        "tema": "Métodos / mates",
    },
    "predecir log raiz": {
        "title": "predecir log raiz",
        "blurb": "Predicción de log y raíz (Newton, comparativas C++/Python, benches).",
        "tema": "Métodos / mates",
    },
    "raiz": {
        "title": "raiz",
        "blurb": "Método de los a para raíz (documento).",
        "tema": "Métodos / mates",
    },
    "riemann": {
        "title": "riemann",
        "blurb": "Notas sobre Riemann / función ζ.",
        "tema": "Métodos / mates",
    },
    "log": {
        "title": "log",
        "blurb": "Logs y experimentos de logaritmo (Python + docx).",
        "tema": "Métodos / mates",
    },
    "COMPARATIVA": {
        "title": "COMPARATIVA — métodos frente a frente",
        "blurb": "Comparativas de cribas, Newton Rápido y raíces (docx + C++).",
        "tema": "Métodos / mates",
    },
    "mates": {
        "title": "mates — matemáticas generales",
        "blurb": "Documentos matemáticos generales VMA (docx/pdf/txt).",
        "tema": "Métodos / mates",
    },
    "Libro1 Números i numeritos": {
        "title": "Libro 1 — Números i numeritos",
        "blurb": "Primer libro del corpus VMA (PDF).",
        "tema": "Libros",
    },
    "Libro2 Números otra vez": {
        "title": "Libro 2 — Números otra vez",
        "blurb": "Segundo libro del corpus VMA (PDF).",
        "tema": "Libros",
    },
    "Libro3 Sigo en mis trece": {
        "title": "Libro 3 — Sigo en mis trece",
        "blurb": "Tercer libro: docx/pdf/txt del corpus.",
        "tema": "Libros",
    },
    "Libro4 Semillas y Energia Verde": {
        "title": "Libro 4 — Semillas y Energía Verde",
        "blurb": "Cuarto libro: semillas, energía verde, ZypyZape / convergencias.",
        "tema": "Libros / energía",
    },
    "Libro5 Factorizacion con 2v+3": {
        "title": "Libro 5 — Factorización con 2v+3",
        "blurb": "Quinto libro: MDC, factorización, scripts Python, gráficas y corpus L5.",
        "tema": "Libros / MDC",
    },
    "Libro6 NewtonRapido en busqueda con oraculo": {
        "title": "Libro 6 — Newton Rápido con oráculo",
        "blurb": "Sexto libro: Newton Rápido y búsqueda con oráculo MEcuation.",
        "tema": "Libros / Newton",
    },
    "libro 4 preparacion": {
        "title": "libro 4 preparacion",
        "blurb": "Borradores y preparación del Libro 4.",
        "tema": "Libros",
    },
    "libro 4 preparacion filesgrans": {
        "title": "libro 4 preparacion filesgrans",
        "blurb": "Archivos grandes de preparación del Libro 4.",
        "tema": "Libros",
    },
    "nuevo libro3": {
        "title": "nuevo libro3",
        "blurb": "Borradores y anexos del Libro 3.",
        "tema": "Libros",
    },
    "pdf de 3 1º libros": {
        "title": "pdf de 3 1º libros",
        "blurb": "PDFs consolidados de los primeros libros (Apiñon, tres…).",
        "tema": "Libros",
    },
    "Archivos del doc Apiñon . docx previo a libro4": {
        "title": "Apiñon — previo a Libro 4",
        "blurb": "Documentos Apiñon y txt previos al Libro 4.",
        "tema": "Libros",
    },
    "libro-metodos-simples": {
        "title": "libro-metodos-simples",
        "blurb": "Notas / libro de métodos simples.",
        "tema": "Libros",
    },
    "PY L5": {
        "title": "PY L5 — Python del Libro 5",
        "blurb": "Scripts Python del Libro 5 / factorización MDC.",
        "tema": "Libros / MDC",
    },
    "files l5": {
        "title": "files l5",
        "blurb": "Auxiliares L5 (ksweep predictivo, tests).",
        "tema": "Libros / MDC",
    },
    "filestot l5": {
        "title": "filestot l5 — total L5",
        "blurb": "Conjunto total de files del Libro 5 (Python, CSV, PNG).",
        "tema": "Libros / MDC",
    },
    "tot l5": {
        "title": "tot l5",
        "blurb": "Consolidado de materiales L5 (espejo de filestot).",
        "tema": "Libros / MDC",
    },
    "txt l5": {
        "title": "txt l5",
        "blurb": "Textos y chats del Libro 5 (incluye corpus Aleatorovix).",
        "tema": "Libros / MDC",
    },
    "DOCX l5": {
        "title": "DOCX l5",
        "blurb": "Documentos Word del Libro 5.",
        "tema": "Libros",
    },
    "FOTOS l5": {
        "title": "FOTOS l5",
        "blurb": "Imágenes, notebooks y capturas del Libro 5.",
        "tema": "Libros",
    },
    "html l5": {
        "title": "html l5",
        "blurb": "HTML y PDF del Libro 5.",
        "tema": "Libros",
    },
    "filesclaude 6-5": {
        "title": "filesclaude 6-5 — exportes Claude",
        "blurb": "Exportes de chat y código (PNG/Python/MD) de sesiones Claude.",
        "tema": "Sesiones IA / archivo",
    },
    "deepseekjun26": {
        "title": "deepseekjun26",
        "blurb": "Sesiones y tablas DeepSeek de junio 2026 (MDC / cribas).",
        "tema": "Sesiones IA / archivo",
    },
    "CLAUDE TXT": {
        "title": "CLAUDE TXT",
        "blurb": "Corpus de chats y PDF asociados a Claude.",
        "tema": "Sesiones IA / archivo",
    },
    "Lee Chat": {
        "title": "Lee Chat",
        "blurb": "Chats exportados (Lee / conversaciones de trabajo).",
        "tema": "Sesiones IA / archivo",
    },
    "lee arbusto": {
        "title": "lee arbusto",
        "blurb": "Material Lee Arbusto: Libro 4 convergencias, zips y textos.",
        "tema": "Sesiones IA / archivo",
    },
    "gemini": {
        "title": "gemini",
        "blurb": "Material de sesiones Gemini (txt/png).",
        "tema": "Sesiones IA / archivo",
    },
    "gemini google comslsshare sls58220b9d68f8": {
        "title": "gemini share export",
        "blurb": "Export HTML/docs de una sesión Gemini compartida (Kuramoto/Zipi-Zape y anexos).",
        "tema": "Sesiones IA / archivo",
    },
    "Kuramoto, Sincronización y Zipi-Zape - Google Gemini_files": {
        "title": "Kuramoto × Zipi-Zape (Gemini files)",
        "blurb": "Assets de la conversación Gemini sobre Kuramoto, sincronización y ZypyZape.",
        "tema": "ZypyZape / sesiones",
    },
    "grok": {
        "title": "grok",
        "blurb": "Material generado o trabajado con Grok (bat, md, txt).",
        "tema": "Sesiones IA / archivo",
    },
    "grokbash": {
        "title": "grokbash",
        "blurb": "Automatizaciones de sesión Grok (jsonl/json de logs).",
        "tema": "Sesiones IA / archivo",
    },
    "copilot 2025": {
        "title": "copilot 2025",
        "blurb": "Sesiones y código con Copilot 2025 (docx/png/pdf).",
        "tema": "Sesiones IA / archivo",
    },
    "New conversation - Grok.html. ojalamia_files (L)": {
        "title": "New conversation Grok (export)",
        "blurb": "Export grande de conversación Grok (html/files).",
        "tema": "Sesiones IA / archivo",
    },
    "gptcomputing": {
        "title": "gptcomputing",
        "blurb": "Experimentos GPT × computación.",
        "tema": "Sesiones IA / archivo",
    },
    "ideas-para-gpt-antipc": {
        "title": "ideas-para-gpt-antipc",
        "blurb": "Brainstorming, roadmaps y specs históricas de AntiPC para GPT.",
        "tema": "AntiPC / archivo",
    },
    "ZYPYZAPE Bateria Cinetica": {
        "title": "ZYPYZAPE Bateria Cinetica",
        "blurb": "ZypyZape: inercia sintética y batería cinética para redes renovables (Python, docs, PNG).",
        "tema": "Energía / ZypyZape",
    },
    "zypyzape quijote ballant": {
        "title": "zypyzape quijote ballant",
        "blurb": "Cruce ZypyZape × Quijote (metáfora del ball de pesos + control).",
        "tema": "Energía / Quijote",
    },
    "zypyzape-contexto": {
        "title": "zypyzape-contexto",
        "blurb": "Contexto y notas alrededor de ZypyZape.",
        "tema": "Energía / ZypyZape",
    },
    "Quijote": {
        "title": "Quijote",
        "blurb": "Quijote: pesos en aspas, hurto gravitatorio, simulaciones y papers.",
        "tema": "Energía / Quijote",
    },
    "Quijotee": {
        "title": "Quijotee",
        "blurb": "Variante / anexos Quijote (más PNG y scripts).",
        "tema": "Energía / Quijote",
    },
    "kilometre;(soles_bateria)": {
        "title": "kilometre (soles_bateria)",
        "blurb": "Kilómetre: máquina de peso rotacional, soles/batería, gráficas y Python.",
        "tema": "Energía / Kilómetre",
    },
    "hurto-gravitatorio": {
        "title": "hurto-gravitatorio",
        "blurb": "Concepto de hurto gravitatorio (metáfora física-numérica del trabajo por revolución).",
        "tema": "Energía",
    },
    "gemelos": {
        "title": "gemelos",
        "blurb": "Gemelos de simulación (primos gemelos y/o gemelos energéticos — revisar subficheros).",
        "tema": "Energía / mates",
    },
    "just run": {
        "title": "just run — listos para ejecutar",
        "blurb": "Paquetes y scripts listos para ejecutar (incluye aleatorovix y gemelos).",
        "tema": "Ejecutables",
    },
    "vma-run": {
        "title": "vma-run — lanzadores",
        "blurb": "Lanzadores y entornos de ejecución VMA (PyInstaller/toc).",
        "tema": "Ejecutables",
    },
    "bin": {
        "title": "bin",
        "blurb": "Binarios o stubs de compilación (p. ej. vma-run.exe).",
        "tema": "Ejecutables",
    },
    "vma": {
        "title": "vma — núcleo VMA",
        "blurb": "Núcleo y materiales VMA (K3 cross-verifier, packs, txt/py).",
        "tema": "VMA núcleo",
    },
    "VMA-Dossier-Elewit": {
        "title": "VMA-Dossier-Elewit",
        "blurb": "Dossier VMA orientado a Elewit / interlocutores industriales.",
        "tema": "VMA / industrial",
    },
    "sale-it": {
        "title": "sale-it — go-to-market",
        "blurb": "Material de venta / go-to-market (scripts, bat, textos).",
        "tema": "VMA / industrial",
    },
    "patentes nacionales": {
        "title": "patentes nacionales",
        "blurb": "Borradores y material de patentes nacionales.",
        "tema": "VMA / legal",
    },
    "seguridad en vehiculos y computacion": {
        "title": "seguridad en vehiculos y computacion",
        "blurb": "Notas de seguridad en vehículos y computación (txt).",
        "tema": "Ingeniería",
    },
    "tecnologia": {
        "title": "tecnologia",
        "blurb": "Notas de tecnología general (docx/pdf).",
        "tema": "Ingeniería",
    },
    "docs": {
        "title": "docs — documentación del monorepo",
        "blurb": "Documentación general: física ZZ+Quijote+Kilómetro, orígenes, integración AntiPC, firma geométrica.",
        "tema": "Docs",
    },
    "biblio": {
        "title": "biblio",
        "blurb": "Bibliografía y referencias (headers C++ / textos).",
        "tema": "Docs",
    },
    "Diapositivas": {
        "title": "Diapositivas",
        "blurb": "Presentaciones y diapositivas del proyecto.",
        "tema": "Docs",
    },
    "Propias": {
        "title": "Propias",
        "blurb": "Presentaciones propias: SaltoMáximo y Siguiente primo (pptx).",
        "tema": "Docs",
    },
    "Filosofía": {
        "title": "Filosofía",
        "blurb": "PDFs de filosofía del autor (marco 33×1 y pensamiento).",
        "tema": "Docs",
    },
    "graficas y explicaciones": {
        "title": "graficas y explicaciones",
        "blurb": "Gráficas MDC/MRAUV/AntiPC y textos explicativos (png/gif/txt).",
        "tema": "Docs / gráficas",
    },
    "codigos gen graficas": {
        "title": "codigos gen graficas",
        "blurb": "Código Python para generar gráficas del corpus.",
        "tema": "Docs / gráficas",
    },
    "fotos": {
        "title": "fotos",
        "blurb": "Fotografías del proyecto (jpg).",
        "tema": "Media",
    },
    "music": {
        "title": "music",
        "blurb": "Audio / música asociada al monorepo.",
        "tema": "Media",
    },
    "py": {
        "title": "py — scripts Python generales",
        "blurb": "Colección de scripts Python y notas md/tex del monorepo.",
        "tema": "Código",
    },
    "pyy": {
        "title": "pyy — más Python",
        "blurb": "Variantes y scripts Python adicionales.",
        "tema": "Código",
    },
    "UNION": {
        "title": "UNION — packs unidos",
        "blurb": "Unión de packs logbench / Docker / benches (cpp, csv, zip).",
        "tema": "Código / benches",
    },
    "conjunto": {
        "title": "conjunto",
        "blurb": "Conjuntos y estructuras (material corto).",
        "tema": "Métodos / mates",
    },
    "diamante": {
        "title": "diamante",
        "blurb": "Proyecto diamante: txt/py/c de exploración numérica / metáfora.",
        "tema": "Experimentos",
    },
    "inicio": {
        "title": "inicio",
        "blurb": "Materiales de punto de entrada / arranque del corpus.",
        "tema": "Docs",
    },
    "archivos-vma": {
        "title": "archivos-vma",
        "blurb": "Archivos varios VMA (miscelánea de proyecto).",
        "tema": "Archivo",
    },
    "desktop-snapshot": {
        "title": "desktop-snapshot",
        "blurb": "Instantánea de archivos clave del Escritorio (HTML Aleatorovix, mapas, guías).",
        "tema": "Archivo",
    },
    "Wolfram": {
        "title": "Wolfram",
        "blurb": "Notas estilo Wolfram / Alpha.",
        "tema": "Sesiones / mates",
    },
    "txt25 5": {
        "title": "txt25 5",
        "blurb": "Textos adicionales línea 5 / 25 (txt/docx/pdf).",
        "tema": "Archivo",
    },
    "junto": {
        "title": "junto — cribas y cotas (docs unidos)",
        "blurb": "Documentos unidos de cribas, cotas y estructuras de primos VMA (docx).",
        "tema": "Métodos / mates",
    },
    "nuevo": {
        "title": "nuevo — cribas/cotas recientes",
        "blurb": "Material reciente: cribas, cotas y estructuras de primos (docx) pendiente de reordenar.",
        "tema": "Métodos / mates",
    },
    "cosas": {
        "title": "cosas — notas sueltas 33×1 / Quijote",
        "blurb": "Notas sueltas: 2a+3, 33×1, Quijote+ZypyZape bilingüe, chats DeepSeek.",
        "tema": "Archivo de trabajo",
    },
    "coses velles": {
        "title": "coses velles — archivo histórico",
        "blurb": "Material antiguo a revisar: Sophie, Criva, docs varios (histórico, no canónico).",
        "tema": "Archivo histórico",
    },
    "otro": {
        "title": "otro — borradores varios",
        "blurb": "Borradores sin clasificar (besons, envíos, páginas VMA).",
        "tema": "Archivo",
    },
    "Nueva carpeta": {
        "title": "Nueva carpeta — borradores libros / primos",
        "blurb": "Borradores de números/primos y PDFs de primeros libros. (No incluye listado de chats personales.)",
        "tema": "Archivo / sin clasificar",
        "privacy": True,
    },
    "Nueva carpeta (3)": {
        "title": "Nueva carpeta (3) — análisis criba Suntaman",
        "blurb": "Análisis de criba Suntaman y librería C (docx).",
        "tema": "Métodos / mates",
    },
    "Nueva carpeta0": {
        "title": "Nueva carpeta0 — criba Suntaman / arquitectura",
        "blurb": "Análisis criba Suntaman y arquitectura (docx).",
        "tema": "Métodos / mates",
    },
    "xd": {
        "title": "xd — Newton / logbench (cajón útil)",
        "blurb": "Cajón con Newton Rápido, benches C++/Python, zips y subcarpeta predecir log raiz.",
        "tema": "Métodos / benches",
    },
    "_extract_vma": {
        "title": "_extract_vma",
        "blurb": "Extracción temporal de packs VMA (no usar como fuente canónica).",
        "tema": "Temporal",
    },
    "tools": {
        "title": "tools — utilidades del monorepo",
        "blurb": "Scripts de mantenimiento: `describe_folders.py` regenera MAPA.md, FOLDERS.md y READMEs por carpeta.",
        "tema": "Docs / tooling",
    },
}

# orden temático para MAPA.md
THEME_ORDER = [
    ("Empezar aquí / 33×1", ["33x1", "33x1-pack", "33x1_CIENCIA", "docs", "VMA_mates_rescat_2026"]),
    ("Aleatorovix + web techamv", ["aleatorovix", "web-aleatorovix", "techamv-aleatorovix-DEPLOY", "techamv-web", "webtechamv", "desktop-snapshot"]),
    ("AntiPC", ["antipc", "antipc-port-c", "antipc2", "cpu-antipc", "ideas-para-gpt-antipc", "sale-it"]),
    ("K3 / cifrado geométrico", ["encriptacionGeometrica", "vma-k3", "hashtool-work", "hashtool-extract", "vma"]),
    ("Métodos y mates", [
        "vma-methods", "Metodo Newton Rápido", "Metodo densidad MRAUV", "mrauv",
        "teoremas", "teoremasgrok", "Goldbach", "Sophie Germain", "sofi",
        "dos_primos", "siguiente_primo", "salto", "discriminant", "densidad",
        "predecir log raiz", "raiz", "riemann", "log", "COMPARATIVA", "mates",
        "conjunto", "junto", "nuevo", "Nueva carpeta (3)", "Nueva carpeta0", "xd",
    ]),
    ("Libros VMA (1–6)", [
        "Libro1 Números i numeritos", "Libro2 Números otra vez", "Libro3 Sigo en mis trece",
        "Libro4 Semillas y Energia Verde", "Libro5 Factorizacion con 2v+3",
        "Libro6 NewtonRapido en busqueda con oraculo",
        "libro 4 preparacion", "libro 4 preparacion filesgrans", "nuevo libro3",
        "pdf de 3 1º libros", "Archivos del doc Apiñon . docx previo a libro4",
        "libro-metodos-simples", "PY L5", "files l5", "filestot l5", "tot l5",
        "txt l5", "DOCX l5", "FOTOS l5", "html l5", "txt25 5",
    ]),
    ("Energía: ZypyZape · Quijote · Kilómetre", [
        "ZYPYZAPE Bateria Cinetica", "zypyzape-contexto", "zypyzape quijote ballant",
        "Quijote", "Quijotee", "kilometre;(soles_bateria)", "hurto-gravitatorio",
        "gemelos", "2026", "1",
    ]),
    ("Ejecutables y lanzadores", ["just run", "vma-run", "bin", "UNION"]),
    ("Sesiones IA / chats de trabajo", [
        "CLAUDE TXT", "Lee Chat", "lee arbusto", "deepseekjun26", "filesclaude 6-5",
        "gemini", "gemini google comslsshare sls58220b9d68f8",
        "Kuramoto, Sincronización y Zipi-Zape - Google Gemini_files",
        "grok", "grokbash", "copilot 2025",
        "New conversation - Grok.html. ojalamia_files (L)", "gptcomputing", "Wolfram",
    ]),
    ("Docs, media, tooling e industrial", [
        "biblio", "Diapositivas", "Propias", "Filosofía", "graficas y explicaciones",
        "codigos gen graficas", "fotos", "music", "inicio", "tecnologia",
        "seguridad en vehiculos y computacion", "patentes nacionales",
        "VMA-Dossier-Elewit", "archivos-vma", "py", "pyy", "tools",
    ]),
    ("Archivo / sin clasificar (histórico)", [
        "cosas", "coses velles", "otro", "Nueva carpeta", "_extract_vma",
    ]),
]

PRIVACY_NAME = re.compile(
    r"(whatsapp|chat de |nif|dni|\+\d{1,3}[\s-]?\d|\b\d{9}\b|contraseña|password|token|secret)",
    re.I,
)
AUTO_MARK = "README generado para navegación del monorepo"
OWN_HINTS = (
    "Adaptive Network",
    "Investigador independent",
    "Qui soc",
    "Fórmula maestra",
    "CONDICIÓ 1",
    "Autor material:",
)


def is_own_readme(text: str) -> bool:
    if AUTO_MARK in text:
        return False
    if len(text) > 1800:
        return True
    return any(h in text for h in OWN_HINTS)


def is_auto_readme(text: str) -> bool:
    return AUTO_MARK in text or text.strip().startswith("# ") and "Monorepo:" in text


def safe_list_files(names: list[str], privacy: bool = False) -> list[str]:
    out = []
    for n in names:
        if n.startswith("~$"):
            continue
        if PRIVACY_NAME.search(n):
            continue
        if privacy and n.lower().endswith((".txt",)) and "chat" in n.lower():
            continue
        out.append(n)
    return out


def inventory(folder: Path) -> dict:
    """Inventario: listado de raíz + conteo total recursivo (aprox.)."""
    files, dirs = [], []
    try:
        for p in folder.iterdir():
            if p.name in {".git", "__pycache__", ".venv", "node_modules"}:
                continue
            if p.is_dir():
                dirs.append(p.name)
            elif p.is_file() and p.name != "README.md":
                files.append(p.name)
    except OSError:
        pass
    files.sort(key=str.lower)
    dirs.sort(key=str.lower)
    exts = Counter()
    total = 0
    try:
        for p in folder.rglob("*"):
            if not p.is_file():
                continue
            parts = set(p.parts)
            if parts & {".git", "__pycache__", ".venv", "node_modules"}:
                continue
            if p.name == "README.md" and p.parent == folder:
                continue
            total += 1
            suf = p.suffix.lower() or "(sin ext)"
            exts[suf] += 1
            if total >= 50000:
                break
    except OSError:
        total = len(files)
        for n in files:
            suf = Path(n).suffix.lower() or "(sin ext)"
            exts[suf] += 1
    return {"files": files, "dirs": dirs, "exts": exts, "total": total}


def describe(name: str, inv: dict) -> dict:
    c = CURATED.get(name)
    if c:
        return {
            "title": c["title"],
            "blurb": c["blurb"],
            "tema": c["tema"],
            "privacy": c.get("privacy", False),
        }
    # heurística por nombre / extensiones
    low = name.lower()
    tema = "Monorepo VMA"
    blurb = f"Carpeta `{name}/` del monorepo. Ver listado de archivos abajo."
    if "libro" in low:
        tema, blurb = "Libros", f"Material de libro / preparación: **{name}**."
    elif "quijote" in low or "zypy" in low or "kilomet" in low:
        tema, blurb = "Energía", f"Material energético / control: **{name}**."
    elif "antipc" in low:
        tema, blurb = "AntiPC", f"Material AntiPC: **{name}**."
    elif any(x in low for x in ("prim", "criba", "newton", "mdc", "gold", "sofi")):
        tema, blurb = "Métodos / mates", f"Material matemático: **{name}**."
    return {"title": name, "blurb": blurb, "tema": tema, "privacy": False}


def render_readme(name: str, meta: dict, inv: dict) -> str:
    privacy = meta.get("privacy", False)
    files = safe_list_files(inv["files"], privacy=privacy)
    dirs = [d for d in inv["dirs"] if d not in {"__pycache__"}]
    exts = inv["exts"]
    n_root = len(inv["files"])
    n_total = inv.get("total", n_root)
    lines = [
        f"# {meta['title']}",
        "",
        meta["blurb"],
        "",
        f"**Tema:** {meta['tema']}",
        "",
        f"**Monorepo:** [{REPO_URL.split('/')[-2]}/{REPO_URL.split('/')[-1]}]({REPO_URL}) · ruta: `{name}/`",
        "",
        "## Contenido",
        "",
        f"- **Archivos totales (aprox., recursivo):** {n_total}",
        f"- **En la raíz de la carpeta:** {n_root}",
        f"- **Subcarpetas (primer nivel):** {len(dirs)}",
        "",
        "### Extensiones (recursivo, top)",
    ]
    if exts:
        for ext, cnt in exts.most_common(12):
            lines.append(f"- `{ext}`: {cnt}")
    else:
        lines.append("_sin archivos listables_")
    lines += ["", "### Subcarpetas"]
    if dirs:
        for d in dirs[:40]:
            lines.append(f"- `{d}`")
        if len(dirs) > 40:
            lines.append(f"- … y {len(dirs) - 40} más")
    else:
        lines.append("_ninguna_")
    lines += ["", "### Archivos en la raíz de esta carpeta"]
    show = files[:40]
    if privacy:
        lines.append("_Listado filtrado (sin chats personales ni datos sensibles)._")
    if show:
        for f in show:
            lines.append(f"- `{f}`")
        if len(files) > 40:
            lines.append(f"- … y {len(files) - 40} más")
    else:
        lines.append("_ninguno listable_ o solo subcarpetas")
    lines += [
        "",
        "## Notas",
        "",
        "- Parte del ecosistema **VMA / 33×1 / AntiPC / Aleatorovix / ZypyZape / K3**.",
        "- Prioridad del monorepo: **33×1** (el **1** = todo el repo civil).",
        "- No subir secretos (claves FTP, tokens, datos personales) a commits futuros.",
        "",
        "## Enlaces relacionados",
        "",
        "- Mapa temático: [MAPA.md](../MAPA.md)",
        "- Índice alfabético: [FOLDERS.md](../FOLDERS.md)",
        "- README raíz: [README.md](../README.md)",
        "- 33×1: [33x1/](../33x1/)",
        "",
        "---",
        f"*{AUTO_MARK}. Editar libremente.*",
        "",
    ]
    return "\n".join(lines)


def md_link(name: str) -> str:
    from urllib.parse import quote

    return f"[{name}]({quote(name)}/README.md)"


def write_mapa(folders: list[str], metas: dict[str, dict], invs: dict[str, dict]) -> None:
    lines = [
        "# Mapa del monorepo `claude`",
        "",
        "**Autor:** Víctor Manzanares Alberola · espiradesombra · VMA",
        "",
        "> **33×1** = cambiar el **1** (TODO este repositorio técnico civil) por **33** años de paz **firmada por los países**.",
        ">",
        "> Empieza por [33x1/00_QUE_ES_33x1.md](33x1/00_QUE_ES_33x1.md) y [33x1/02_USO_CIVIL.txt](33x1/02_USO_CIVIL.txt).",
        "",
        "Este mapa agrupa las carpetas de **primer nivel** por tema. Cada carpeta tiene un `README.md` con su contenido.",
        "Índice alfabético completo: [FOLDERS.md](FOLDERS.md).",
        "",
    ]
    placed: set[str] = set()
    for theme, names in THEME_ORDER:
        present = [n for n in names if n in metas]
        if not present:
            continue
        lines.append(f"## {theme}")
        lines.append("")
        lines.append("| Carpeta | Archivos | Qué contiene |")
        lines.append("|---------|----------|--------------|")
        for name in present:
            inv = invs[name]
            blurb = metas[name]["blurb"].replace("|", "\\|")
            if len(blurb) > 140:
                blurb = blurb[:137] + "…"
            count = inv.get("total", len(inv["files"]))
            lines.append(f"| {md_link(name)} | {count}+ | {blurb} |")
            placed.add(name)
        lines.append("")
    rest = sorted([n for n in folders if n not in placed], key=str.lower)
    if rest:
        lines.append("## Otras carpetas")
        lines.append("")
        lines.append("| Carpeta | Archivos | Qué contiene |")
        lines.append("|---------|----------|--------------|")
        for name in rest:
            inv = invs[name]
            blurb = metas[name]["blurb"].replace("|", "\\|")
            if len(blurb) > 140:
                blurb = blurb[:137] + "…"
            count = inv.get("total", len(inv["files"]))
            lines.append(f"| {md_link(name)} | {count}+ | {blurb} |")
        lines.append("")
    lines += [
        "---",
        "",
        "## Convención de commits (regla A)",
        "",
        "**Un commit = una carpeta.** El mensaje nombra la carpeta y describe contenido o cambio.",
        "Detalle: [COMMITS.md](COMMITS.md).",
        "",
        "```",
        "docs(VMA_mates_rescat_2026): pack mates — cribas fix, MDC, XFI N=3",
        "feat(web-aleatorovix): motor JS sin Math.random + UI del ramo",
        "docs(raiz): MAPA + FOLDERS",
        "```",
        "",
        "Regenerar este mapa:",
        "",
        "```bat",
        "python tools\\describe_folders.py",
        "```",
        "",
    ]
    (ROOT / "MAPA.md").write_text("\n".join(lines), encoding="utf-8")


def write_folders(folders: list[str], metas: dict[str, dict], invs: dict[str, dict]) -> None:
    lines = [
        "# Índice de carpetas — monorepo `claude`",
        "",
        "Autor: Víctor Manzanares Alberola · espiradesombra · VMA / 33×1",
        "",
        "Cada carpeta de primer nivel tiene un `README.md` que describe **la carpeta y su contenido**.",
        "Mapa temático (recomendado): **[MAPA.md](MAPA.md)**.",
        "",
        f"**Carpetas documentadas:** {len(folders)}",
        "",
        "| Carpeta | Archivos | Tema | Resumen |",
        "|---------|----------|------|---------|",
    ]
    for name in sorted(folders, key=str.lower):
        inv = invs[name]
        m = metas[name]
        blurb = m["blurb"].replace("|", "\\|")
        if len(blurb) > 100:
            blurb = blurb[:97] + "…"
        tema = m["tema"].replace("|", "\\|")
        count = inv.get("total", len(inv["files"]))
        lines.append(f"| {md_link(name)} | {count}+ | {tema} | {blurb} |")
    lines += [
        "",
        "---",
        "",
        "Generado con `python tools/describe_folders.py`.",
        "",
    ]
    (ROOT / "FOLDERS.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="Reescribe también READMEs propios/manuales")
    ap.add_argument("--force-weak", action="store_true", help="Reescribe solo auto o genéricos cortos")
    ap.add_argument("--only", type=str, default="", help="Solo estas carpetas (coma-separadas)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    only = {x.strip() for x in args.only.split(",") if x.strip()}

    folders = sorted(
        [p.name for p in ROOT.iterdir() if p.is_dir() and not p.name.startswith(".git")],
        key=str.lower,
    )
    metas: dict[str, dict] = {}
    invs: dict[str, dict] = {}
    written = 0
    skipped = 0

    for name in folders:
        path = ROOT / name
        inv = inventory(path)
        invs[name] = inv
        meta = describe(name, inv)
        metas[name] = meta

        if only and name not in only:
            continue

        readme_path = path / "README.md"
        existing = ""
        if readme_path.exists():
            existing = readme_path.read_text(encoding="utf-8", errors="replace")

        should = False
        if args.force:
            should = True
        elif not existing:
            should = True
        elif args.force_weak and (is_auto_readme(existing) or len(existing) < 900):
            # no pisar propios largos
            if not is_own_readme(existing):
                should = True
        elif is_auto_readme(existing):
            should = True
        elif name in CURATED and is_auto_readme(existing):
            should = True

        # siempre reescribir auto + curated cuando force-weak o default auto
        if not should and is_auto_readme(existing) and name in CURATED:
            should = True

        if not should:
            skipped += 1
            continue

        body = render_readme(name, meta, inv)
        if args.dry_run:
            print(f"WOULD WRITE {name}/README.md")
        else:
            readme_path.write_text(body, encoding="utf-8")
            print(f"WRITE {name}/README.md")
        written += 1

    if not args.dry_run:
        write_mapa(folders, metas, invs)
        write_folders(folders, metas, invs)
        print(f"MAPA.md + FOLDERS.md actualizados ({len(folders)} carpetas)")
    print(f"readmes written={written} skipped_own_or_filtered={skipped}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
