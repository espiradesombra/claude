// El hash cambia con el tiempo basado en su propio valor
function hashAutoModificable(input) {
    let hash = calcularHash(input);
    // El hash se usa para modificar sus propios parámetros
    const nuevoFactor = 1 + (hash % 10) / 10;
    return calcularHashConFactor(input, nuevoFactor);
}