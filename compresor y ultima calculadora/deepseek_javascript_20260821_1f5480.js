// Los punteros se mueven con reglas fijas
while (i <= j) {
    v *= factor1;  // Acelera
    i += punteroA;
    j -= punteroB;
    if (i > j) { choque(); }  // ¡Caos!
}