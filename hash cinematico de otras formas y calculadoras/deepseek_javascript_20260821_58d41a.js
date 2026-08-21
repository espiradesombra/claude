// Sistema completo de seguridad por inferencias de ondas
class SistemaSeguridadOndas {
    constructor() {
        this.historial = [];
        this.parametros = null;
        this.hash = null;
    }
    
    // Simular captura de ondas
    capturarOndas() {
        // Simular dos frecuencias
        const csi2_4 = capturarCSI(2.4e9);
        const csi5 = capturarCSI(5e9);
        
        // Extraer parámetros físicos
        const datosFisicos = extraerParametros(csi2_4, csi5);
        
        // Guardar historial
        this.historial.push({
            timestamp: Date.now(),
            datos: datosFisicos,
            csi2_4,
            csi5
        });
        
        return datosFisicos;
    }
    
    // Generar hash de seguridad
    generarHash() {
        // Capturar ondas
        const datosFisicos = this.capturarOndas();
        
        // Generar parámetros desde las ondas
        this.parametros = generarParametrosDesdeOndas(datosFisicos);
        
        // Crear string para hash
        const hashData = [
            `distancia=${datosFisicos.distancia}`,
            `velocidad=${datosFisicos.velocidad}`,
            `doppler=${datosFisicos.doppler}`,
            `multicamino=${JSON.stringify(datosFisicos.multicamino)}`,
            `timestamp=${datosFisicos.timestamp}`,
            `parametros=${JSON.stringify(this.parametros)}`
        ];
        
        const fullStr = hashData.join('|');
        const bigInt = stringToBigInt(fullStr);
        const mod = 2n ** BigInt(64);
        this.hash = bigInt % mod;
        
        return {
            hash: bigIntToHex(this.hash, 64),
            parametros: this.parametros,
            datosFisicos: datosFisicos,
            historial: this.historial
        };
    }
    
    // Verificar si el entorno ha cambiado
    verificarEntorno() {
        if (this.historial.length < 2) {
            return { estable: false, mensaje: 'Se necesitan más mediciones' };
        }
        
        const ultimo = this.historial[this.historial.length - 1];
        const anterior = this.historial[this.historial.length - 2];
        
        const cambioDistancia = Math.abs(ultimo.datos.distancia - anterior.datos.distancia);
        const cambioVelocidad = Math.abs(ultimo.datos.velocidad - anterior.datos.velocidad);
        
        return {
            estable: cambioDistancia < 0.5 && cambioVelocidad < 0.3,
            cambioDistancia,
            cambioVelocidad,
            mensaje: cambioDistancia < 0.5 && cambioVelocidad < 0.3 ? 
                '✅ Entorno estable' : '⚠️ Cambios detectados en el entorno'
        };
    }
}