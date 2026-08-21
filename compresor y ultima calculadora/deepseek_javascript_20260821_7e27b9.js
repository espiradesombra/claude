// Reflectividad como función del tiempo
function reflectividad(t, params) {
    const { R0, k, intensidad, angulo } = params;
    // R(t) = R0 * e^(-k * intensidad * t * cos(angulo))
    return R0 * Math.exp(-k * intensidad * t * Math.cos(angulo));
}

// El hash es el estado del espejo en el momento T
function hashEspejo(t, params) {
    const R = reflectividad(t, params);
    // Usar la reflectividad para generar el hash
    return generarHashDesdeReflectividad(R, t, params);
}