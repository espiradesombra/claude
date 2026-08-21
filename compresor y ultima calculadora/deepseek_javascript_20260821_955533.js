function hashUniversal(input, contexto) {
    // Combinar TODAS las fuentes de entropía
    const hash1 = hashWiFi(input, contexto.csi);
    const hash2 = hashAstronomico(input, contexto.planetas);
    const hash3 = hashVoz(input, contexto.audio);
    const hash4 = hashCardiaco(input, contexto.corazon);
    const hash5 = hashOndas(input, contexto.ondasCerebrales);
    
    // Combinar con XOR y rotaciones
    let final = hash1 ^ hash2 ^ hash3 ^ hash4 ^ hash5;
    final = rotar(final, 13);
    final ^= (hash1 << 32) | (hash2 >> 32);
    
    return final;
}