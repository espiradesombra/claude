# Solo los parámetros más importantes (1, 2, 4)
parametros_criticos = {
    'factor1': [1.1, 1.2, 1.3, 1.4, 1.5],      # Sobrepaso
    'factor2': [3.0, 3.2, 3.4, 3.6, 3.8, 4.0],  # Multiplicativo  
    'factor4': [1.5, 1.6, 1.7, 1.8, 1.9, 2.0]   # Exponencial
}

# Mantener fijos los menos críticos
F3_FIJO = -0.85
F5_FIJO = -0.85  
F6_FIJO = 1.8