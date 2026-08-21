// Parámetros de hash desde órbitas
const params = {
    factor1: 1 + (periodo_orbital % 10) / 5,
    factor0: 0.1 + (excentricidad % 10) / 20,
    pasadas: Math.floor(velocidad_escape / 1000)
};