// Esto verá:
const vFlotacion = Math.sqrt(2 * gravedad * altura) * 
    (densidadFluido - densidadObjeto) / densidadFluido * factorFlotacion;

// Y luego:
if (resetVelocidad) {
    vA = vFlotacion;
} else {
    vA += vFlotacion * 0.1;
}

// Y luego C:
const pasoC = vC * 0.2 * direccionC * (1 + Math.sin(t * 0.05) * 0.05);