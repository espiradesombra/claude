// 1. Dividir datos en bloques
const bloques = dividirEnBloques(datos, 16); // 16 bits cada bloque

// 2. Cada bloque → hash reversible
bloques.forEach(bloque => {
    const hash = generarHashReversible(bloque, params);
    // hash = a3f1c78f2b5e4d9c
    // params = { R0: 0.95, k: 0.01, ... }
    guardar(hash, params);
});

// 3. Guardar solo los hashes
archivo_comprimido = {
    hashes: [hash1, hash2, hash3, ...],
    params: [params1, params2, params3, ...]
};

// 4. Descomprimir
archivo_comprimido.hashes.forEach((hash, i) => {
    const params = archivo_comprimido.params[i];
    const bloque = reconstruirBloque(hash, params);
    // bloque = "hola"
});