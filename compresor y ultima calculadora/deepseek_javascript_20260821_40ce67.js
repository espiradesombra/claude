// Usar datos de EEG (ondas cerebrales)
function hashOnirico(input, ondasCerebrales) {
    const alpha = analizarAlpha(ondasCerebrales);
    const beta = analizarBeta(ondasCerebrales);
    const theta = analizarTheta(ondasCerebrales);
    return hashConCerebro(input, { alpha, beta, theta });
}