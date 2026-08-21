// Usar datos de sensores de gas/olor
function hashOlfativo(input, sensores) {
    const compuestos = analizarCompuestos(sensores);
    const concentracion = analizarConcentracion(sensores);
    return hashConOlor(input, { compuestos, concentracion });
}