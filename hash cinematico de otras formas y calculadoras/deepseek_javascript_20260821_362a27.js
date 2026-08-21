// Simulación de captura de CSI para dos frecuencias
function capturarCSI(frecuencia) {
    // Simular parámetros físicos reales
    const distancia = 1 + Math.random() * 5; // 1-5 metros
    const velocidad = Math.random() * 2; // 0-2 m/s
    const angulo = Math.random() * 360; // 0-360 grados
    const multicamino = Array.from({length: 5}, () => ({
        distancia: Math.random() * 10,
        atenuacion: 0.1 + Math.random() * 0.9,
        fase: Math.random() * 2 * Math.PI
    }));
    
    // Calcular fase
    const fase = 2 * Math.PI * frecuencia * distancia / 3e8;
    
    // Calcular Doppler
    const doppler = -2 * velocidad * frecuencia / 3e8 * Math.cos(angulo);
    
    return {
        frecuencia,
        distancia,
        velocidad,
        angulo,
        fase,
        doppler,
        multicamino,
        timestamp: Date.now()
    };
}