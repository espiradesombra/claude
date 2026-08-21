// Usar movimiento del mouse como entropía
function hashBailarin(input, movimientos) {
    const velocidad = calcularVelocidad(movimientos);
    const aceleracion = calcularAceleracion(movimientos);
    const patron = analizarPatron(movimientos);
    return hashConMovimiento(input, { velocidad, aceleracion, patron });
}