"""Métodos de los Libros 1–6 VMA — puente al repo y AntiPC."""

from .bridge import (
    LIBRO_REGISTRY,
    format_libro_info,
    format_metodos_table,
    list_libros,
    list_metodos,
    resolve_metodo,
    run_metodo,
)

__all__ = [
    "LIBRO_REGISTRY",
    "format_libro_info",
    "format_metodos_table",
    "list_libros",
    "list_metodos",
    "resolve_metodo",
    "run_metodo",
]