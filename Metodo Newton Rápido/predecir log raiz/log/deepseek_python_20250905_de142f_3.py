import numpy as np
from itertools import product

# Rangos de prueba para cada parámetro
rangos = {
    'factor1': np.arange(1.1, 1.6, 0.1),      # 1.1, 1.2, 1.3, 1.4, 1.5
    'factor2': np.arange(3.0, 4.1, 0.2),      # 3.0, 3.2, 3.4, 3.6, 3.8, 4.0
    'factor3': np.arange(-1.0, -0.5, 0.1),    # -1.0, -0.9, -0.8, -0.7, -0.6
    'factor4': np.arange(1.5, 2.1, 0.1),      # 1.5, 1.6, 1.7, 1.8, 1.9, 2.0
    'factor5': np.arange(-1.0, -0.5, 0.1),    # -1.0, -0.9, -0.8, -0.7, -0.6  
    'factor6': np.arange(1.5, 2.1, 0.1)       # 1.5, 1.6, 1.7, 1.8, 1.9, 2.0
}

# Casos de prueba diversos
casos_prueba = [
    (100, 10), (8, 2), (27, 3), (1024, 2), 
    (10000, 10), (0.001, 10), (math.pi, 10),
    (123456, 10), (0.123, 10)
]