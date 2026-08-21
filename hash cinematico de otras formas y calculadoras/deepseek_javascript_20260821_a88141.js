// Generación del hash (irreversible sin parámetros)
function generarHash(input, params) {
    // 1. Simulación con parámetros
    const resultado = simularOceanico(input, params);
    
    // 2. Convertir a string científico
    const dataStr = JSON.stringify(resultado);
    
    // 3. Hash final (mod 2^N)
    const bigInt = stringToBigInt(dataStr);
    return bigInt % 2n ** BigInt(params.bits);
}

// Reversión (solo con parámetros)
function revertirHash(hash, params) {
    // 1. Usar los parámetros para reconstruir
    // 2. Buscar el input que produce ese hash
    // 3. (Búsqueda determinista gracias a los parámetros)
    return inputEncontrado;
}