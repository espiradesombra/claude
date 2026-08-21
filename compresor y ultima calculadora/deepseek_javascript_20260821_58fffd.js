// Usar datos de ritmo cardíaco
function hashCardiaco(input, ritmoCardiaco) {
    const variabilidad = analizarVariabilidad(ritmoCardiaco);
    const patron = analizarPatronCardiaco(ritmoCardiaco);
    return hashConCorazon(input, { variabilidad, patron });
}