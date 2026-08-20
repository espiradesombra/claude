def approx_e_custom(skip_modulos, n_terms):
    """
    Aproximación de e por suma de inversos de enteros, descartando múltiplos
    de los módulos en skip_modulos, y sumando 2 al final.
    skip_modulos: lista de enteros cuyos múltiplos se excluyen.
    n_terms: número de términos a considerar en la suma.
    """
    acum = 0.0
    count = 0
    i = 1
    while count < n_terms:
        # Comprobar si i es múltiplo de algún módulo en skip_modulos
        if any(i % m == 0 for m in skip_modulos):
            i += 1
            continue
        acum += 1.0 / i
        count += 1
        i += 1
    return 2.0 + acum