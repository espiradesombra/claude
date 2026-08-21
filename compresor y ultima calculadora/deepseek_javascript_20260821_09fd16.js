// 1. Datos originales: "hola mundo como estas"
// 2. Dividir en bloques: ["hola", "mun", "do ", "com", "o e", "sta", "s"]
// 3. Generar hashes con CLAVE MAESTRA
// 4. Guardar ORDEN: [0:a3, 1:f1, 2:c7, 3:b9, 4:d4, 5:e2, 6:8f]
// 5. Comprimido: "a3f1c7b9d4e28f" (8 bytes)
// 6. ¡Original 28 bytes → Comprimido 8 bytes!
// 7. Descompresión: usa el ORDEN + CLAVE MAESTRA para reconstruir