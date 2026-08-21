// Usar análisis de voz como parámetro
function hashVoz(input, audioData) {
    const frecuencia = analizarFrecuencia(audioData);
    const tono = analizarTono(audioData);
    const ritmo = analizarRitmo(audioData);
    return hashConParametrosAudio(input, { frecuencia, tono, ritmo });
}