// Usar datos de temperatura y humedad
function hashComestible(input, temp, humedad) {
    const indiceCalor = calcularIndiceCalor(temp, humedad);
    const puntoRocio = calcularPuntoRocio(temp, humedad);
    return hashConClima(input, { indiceCalor, puntoRocio });
}