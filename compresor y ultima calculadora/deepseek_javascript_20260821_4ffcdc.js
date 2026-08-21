// Cada fotón es una medición
function medirFoton(estadoEspejo) {
    const probabilidadReflejo = estadoEspejo.R;
    return Math.random() < probabilidadReflejo ? 'reflejado' : 'absorbido';
}

// El hash es el colapso de la función de onda
function colapsarHash(estadoEspejo) {
    // El estado cuántico se "paraliza"
    return generarHash(estadoEspejo);
}