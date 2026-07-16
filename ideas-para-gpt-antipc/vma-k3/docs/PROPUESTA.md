# Propuesta VMA — repo unificado

## Diagnóstico de lo que tenías

| Fuente | Problema | Acción |
|--------|----------|--------|
| `hashcode.zip` / `k3hash.zip` | Duplicados + `has.txt`/`code.txt` enormes mezclando chat y código | Separar motor limpio en `vma/` |
| `code.txt` | Motor industrial C completo pero monolítico | Fase 2: `c/k3_audit/` |
| `has.txt` | Múltiples versiones Python/C/Java repetidas | Una sola `cross_verifier.py` canónica |
| Chats Gemini/Claude | Ideas buenas, dispersas | `docs/CRONOLOGIA-CODIGOS.md` + `glossary.py` |

## Arquitectura propuesta

```
vma/
├── vma/k3/          ← Python canónico (auditoría industrial)
├── c/k3hash/        ← Librería ligera streaming (ya existente)
├── c/k3_audit/      ← [FASE 2] Motor C con firma + stride + marcas
├── examples/        ← Modos Usual / Propio A / Propio B
├── tests/           ← Regresión determinista
└── docs/            ← Especificación y glosario
```

## Por qué empezar por Python

1. El `K3CrossVerifier` del chat ya está completo y probado conceptualmente.
2. Permite validar los 3 modos sin compilar C en Windows.
3. El C de `code.txt` se porta cuando los tests Python fijen hashes de referencia.

## Herramientas/skills que crear después

| Skill | Para qué |
|-------|----------|
| `k3-toolkit` | Compilar, testear, parametrizar desde Grok |
| `zip-evaluator` | Revisar instalables que subas |
| `conceptual-map-of-chat` | Mapas ASCII de sesiones largas como la tuya |

## Prioridad inmediata (ya hecho en v0.1)

- [x] `K3CrossVerifier` Python
- [x] CLI `vma-k3`
- [x] Tests unitarios
- [x] Ejemplos 3 modos
- [x] Copia `c/k3hash` del zip
- [ ] Hashes golden cross Python↔C
- [ ] `k3_audit.c` portado desde `code.txt`
- [ ] pip wheel con binario nativo