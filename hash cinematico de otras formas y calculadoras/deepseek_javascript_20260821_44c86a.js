// Extraer parámetros de las dos frecuencias
function extraerParametros(csi2_4, csi5) {
    // Teorema Chino del Resto para resolver ambigüedad de distancia
    function resolverDistanciaCRT(fase1, fase2, freq1, freq2) {
        // Resolver congruencias: d = (fase_i * c) / (2 * pi * freq_i)
        const d1 = (fase1 * 3e8) / (2 * Math.PI * freq1);
        const d2 = (fase2 * 3e8) / (2 * Math.PI * freq2);
        
        // Encontrar la distancia real usando CRT
        // En la práctica, se busca la solución que minimiza el error
        const posibles = [];
        for (let k = -10; k <= 10; k++) {
            const d = d1 + k * 3e8 / (2 * Math.PI * (freq2 - freq1));
            posibles.push(d);
        }
        
        // Elegir la distancia positiva más cercana
        const distancias = posibles.filter(d => d > 0);
        return distancias.reduce((a, b) => Math.abs(a - d2) < Math.abs(b - d2) ? a : b);
    }
    
    // Resolver distancia
    const distancia = resolverDistanciaCRT(
        csi2_4.fase,
        csi5.fase,
        csi2_4.frecuencia,
        csi5.frecuencia
    );
    
    // Calcular velocidad relativa
    const velocidad = Math.abs(csi2_4.velocidad - csi5.velocidad);
    
    // Calcular Doppler promedio
    const doppler = (csi2_4.doppler + csi5.doppler) / 2;
    
    // Obtener multicamino
    const multicamino = csi2_4.multicamino.map((m, i) => ({
        ...m,
        distancia: m.distancia,
        atenuacion: m.atenuacion * (1 + Math.random() * 0.1)
    }));
    
    return {
        distancia,
        velocidad,
        doppler,
        multicamino,
        timestamp: Date.now()
    };
}