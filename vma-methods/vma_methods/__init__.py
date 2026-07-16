"""
vma-methods — Teoría de números VMA.

  import vma_methods as vma
  vma.cribas.CribaModular6k(5000).run()
  vma.criva.compare_criva_vs_pnt([100, 1000, 10000])
  vma.newton.newton_rapido(121, b=10)
"""

from . import cribas, criva, newton, classic

__version__ = "0.1.0"
__all__ = ["cribas", "criva", "newton", "classic", "__version__"]