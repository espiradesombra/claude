// Generar parámetros para hash cinemático desde datos físicos
function generarParametrosDesdeOndas(datosFisicos) {
    const { distancia, velocidad, doppler, multicamino } = datosFisicos;
    
    // Usar los valores físicos para parametrizar el hash
    const params = {
        // Para hash cinemático
        factor1: 1 + (distancia % 10) / 5,
        factor0: 0.1 + (velocidad % 10) / 20,
        pasadas: Math.max(1, Math.floor(distancia * 2)),
        
        // Para hash nuclear
        lambda1: 0.1 + (doppler % 5) / 5,
        lambda0: 0.05 + (velocidad % 5) / 10,
        lambdaC: 0.5 + (distancia % 10) / 5,
        
        // Para hash oceánico
        gravedad: 1 + (distancia % 20) / 2,
        altura: 10 + (Math.abs(velocidad) % 100) * 10,
        densidadFluido: 0.998 + (Math.abs(doppler) % 100) / 1000,
        
        // Para 3 punteros
        factorA: 1 + (velocidad % 10) / 5,
        factorB: 0.1 + (doppler % 10) / 20,
        factorC: 0.5 + (distancia % 10) / 5,
        vCInicial: 0.1 + (Math.abs(velocidad) % 100) / 100,
        
        // Metadatos de las ondas
        frecuencia2_4: 2.4e9,
        frecuencia5: 5e9,
        distancia: distancia,
        velocidad: velocidad,
        doppler: doppler,
        multicamino: multicamino.map(m => ({
            distancia: m.distancia,
            atenuacion: m.atenuacion
        }))
    };
    
    return params;
}