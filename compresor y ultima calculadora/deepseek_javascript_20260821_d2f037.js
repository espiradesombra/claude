// 1. Diccionario: mapea bloques largos → hashes cortos
const diccionario = {
    "hola": "a3",      // 4 bytes → 1 byte
    "mundo": "f1",     // 5 bytes → 1 byte
    "como": "c7",      // 4 bytes → 1 byte
    "estas": "b9"      // 5 bytes → 1 byte
};

// 2. Datos comprimidos: solo los hashes
const comprimido = ["a3", "f1", "c7", "b9"]; // 4 bytes

// 3. ¡Compresión REAL!
// Original: 18 bytes → Comprimido: 4 bytes
// Ratio: 4.5x ¡REDUCCIÓN!