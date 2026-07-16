"""MDC / Regla Mecánica — port desde HASHTOOLCODE cpp y dll/biblio."""

from .regla_calculo import ReglaCalculo
from .analisi_modular import fase_modular, patron_bits, curvatura_modular
from .factoritzacio import k_sweep_mdc, factorizar_mdc_toy

__all__ = [
    "ReglaCalculo",
    "fase_modular",
    "patron_bits",
    "curvatura_modular",
    "k_sweep_mdc",
    "factorizar_mdc_toy",
]