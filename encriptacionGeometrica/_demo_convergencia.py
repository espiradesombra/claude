"""Demo segura: convergencia_binaria (sin borrar archivos)."""
from decimal import Decimal, getcontext

getcontext().prec = 100


def obtener_perimetro_figura(tipo, puntos_origen, escala_tales):
    base_perimetro = Decimal(puntos_origen) * Decimal("1.5")
    if tipo.lower() == "equilatero":
        lado = (base_perimetro / Decimal("3")) * escala_tales
        return [lado, lado, lado]
    if tipo.lower() == "isosceles":
        lado_a = (base_perimetro / Decimal("4")) * escala_tales
        lado_b = (base_perimetro / Decimal("2")) * escala_tales
        return [lado_a, lado_a, lado_b]
    l1 = (base_perimetro * Decimal("0.25")) * escala_tales
    l2 = (base_perimetro * Decimal("0.35")) * escala_tales
    l3 = (base_perimetro * Decimal("0.40")) * escala_tales
    return [l1, l2, l3]


def encriptar_convergencia(bits, tales, figuras, puntos):
    perimetro_acumulado = Decimal("0.0")
    idx_bit = 0
    iteracion = 0
    n_tales = len(tales)
    n_figuras = len(figuras)

    print("\n--- INICIANDO PROCESO DE ENCRIPTACION GEOMETRICA ---")

    while idx_bit < len(bits):
        par_bits = bits[idx_bit : idx_bit + 2]
        if len(par_bits) == 1:
            par_bits += "0"

        escala = Decimal(str(tales[iteracion % n_tales]))
        figura_tipo = figuras[iteracion % n_figuras]
        puntos_fig = puntos[iteracion % len(puntos)]
        lados = obtener_perimetro_figura(figura_tipo, puntos_fig, escala)
        q = Decimal("1.0000")
        det = (
            f"Iter {iteracion:02d} | Bits: {par_bits} | Tales: {escala} | "
            f"{figura_tipo} ({puntos_fig}pts)"
        )

        if par_bits == "10":
            perimetro_acumulado += lados[0]
            print(f"{det} => APORTA Lado 1 ({lados[0].quantize(q)})")
            idx_bit += 1
        elif par_bits == "11":
            perimetro_acumulado += lados[1]
            print(f"{det} => APORTA Lado 2 ({lados[1].quantize(q)})")
            idx_bit += 1
        elif par_bits == "00":
            perimetro_acumulado += lados[0] + lados[1]
            print(f"{det} => APORTA L1+L2 ({(lados[0]+lados[1]).quantize(q)})")
            idx_bit += 1
        elif par_bits == "01":
            print(f"{det} => APORTA: Ninguno (Salto Doble)")
            idx_bit += 2

        iteracion += 1

    return perimetro_acumulado


if __name__ == "__main__":
    tales = [3, 5, 8, 13, 21]
    figuras = ["equilatero", "isosceles", "escaleno"]
    puntos = [6, 12, 18]
    cadena = "0101101011000010"

    print("Datos de configuracion:")
    print(f"  Tales   : {tales}")
    print(f"  Figuras : {figuras}  puntos {puntos}")
    print(f"  Binario : {cadena}  ({len(cadena)} bits)")

    resultado = encriptar_convergencia(cadena, tales, figuras, puntos)

    print("\n" + "=" * 60)
    print("PERIMETRO DE COLAPSO (encriptado):")
    print(resultado)
    print("=" * 60)