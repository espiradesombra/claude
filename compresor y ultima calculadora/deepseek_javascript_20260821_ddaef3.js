// El hash depende de la posición de los planetas en el futuro/pasado
function hashTemporal(input, offset) {
    const fecha = new Date();
    fecha.setHours(fecha.getHours() + offset);
    const posicionPlanetas = calcularPosiciones(fecha);
    return hashConParametrosAstronomicos(input, posicionPlanetas);
}