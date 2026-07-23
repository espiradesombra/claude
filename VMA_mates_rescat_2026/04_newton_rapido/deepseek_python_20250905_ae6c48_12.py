# Con los mejores parámetros críticos, optimizar los demás
MEJOR_F1, MEJOR_F2, MEJOR_F4 = top_5[0][1], top_5[0][2], top_5[0][3]

# Ahora optimizar factores 3, 5, 6
parametros_secundarios = {
    'factor3': [-0.9, -0.85, -0.8, -0.75],
    'factor5': [-0.9, -0.85, -0.8, -0.75], 
    'factor6': [1.7, 1.8, 1.9, 2.0]
}