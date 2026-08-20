import os
import sys
import secrets
from fractions import Fraction

# ------------------------------------------------------------
# 1. Secuencias y proporciones
# ------------------------------------------------------------

def secuencia_thales(base):
    """Genera secuencia de Thales: base, base+2, base+3, base+5, base+8"""
    return [base, base+2, base+3, base+5, base+8]

# Proporciones de lados para cada tipo de triángulo
RATIOS = {
    1: [Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)],          # Equilátero
    2: [Fraction(1, 4), Fraction(1, 4), Fraction(1, 2)],          # Isósceles
    3: [Fraction(1, 4), Fraction(7, 20), Fraction(2, 5)]          # Escaleno
}

# ------------------------------------------------------------
# 2. Factor geométrico acumulado (según Libro 4)
# ------------------------------------------------------------

def factor_geometrico(base, tipos_triangulo, num_ciclos, V):
    """
    Acumula lados de triángulos escalados por Thales.
    Obtiene a = suma_total % 1  (factor fraccional entre 0 y 1)
    Factor final = (1 + a) * V  (con V como secreto)
    Devuelve el numerador del factor como entero.
    """
    S = secuencia_thales(base)
    a = Fraction(0)

    for i in range(num_ciclos):
        tipo = tipos_triangulo[i % len(tipos_triangulo)]
        escala = S[i % len(S)]
        for r in RATIOS[tipo]:
            a += r * escala   # acumula lados o relaciones

    a = a % 1                # nos quedamos con la parte fraccionaria
    factor_frac = 1 + a      # (1 ± a)
    F = factor_frac.numerator  # tomamos el numerador

    F_total = F * V          # añadimos el factor secreto (1 ± v)
    return F_total

# ------------------------------------------------------------
# 3. Encriptar / desencriptar
# ------------------------------------------------------------

def encriptar_archivo(ruta_entrada, ruta_salida, base, tipos_triangulo, num_ciclos):
    # Leer archivo original
    with open(ruta_entrada, 'rb') as f:
        datos = f.read()

    original_length = len(datos)
    M = int.from_bytes(datos, 'big')          # segmento como entero

    # Generar V aleatorio de 128 bits (factor secreto)
    V = secrets.randbits(128) + 1            # evitamos V=0

    # Calcular factor geométrico total
    F_total = factor_geometrico(base, tipos_triangulo, num_ciclos, V)

    # Cifrar: multiplicar el segmento por el factor
    C = M * F_total

    # Convertir a bytes
    cifrado_bytes = C.to_bytes((C.bit_length() + 7) // 8, 'big') if C > 0 else b'\x00'

    # Guardar archivo con metadatos
    # Formato:
    # magic(4) | version(1) | original_length(8) | num_ciclos(4) | base(8) | V(16)
    # n_tipos(1) | tipos(n_tipos) | cifrado_bytes
    magic = b'GEO1'
    version = 1
    header = magic + bytes([version])
    header += original_length.to_bytes(8, 'big')
    header += num_ciclos.to_bytes(4, 'big')
    header += base.to_bytes(8, 'big')
    header += V.to_bytes(16, 'big')
    header += bytes([len(tipos_triangulo)])
    for t in tipos_triangulo:
        header += bytes([t])
    # Añadir los datos cifrados
    with open(ruta_salida, 'wb') as f:
        f.write(header + cifrado_bytes)

    print(f"✅ Archivo encriptado guardado como: {ruta_salida}")
    print(f"   Tamaño original: {original_length} bytes")
    print(f"   Tamaño cifrado total: {len(header) + len(cifrado_bytes)} bytes")
    print(f"   Factor geométrico F_total: {F_total}")

def desencriptar_archivo(ruta_entrada, ruta_salida):
    with open(ruta_entrada, 'rb') as f:
        data = f.read()

    # Parsear encabezado
    magic = data[0:4]
    if magic != b'GEO1':
        raise ValueError("Formato no reconocido")

    pos = 4
    version = data[pos]; pos += 1
    if version != 1:
        raise ValueError("Versión no soportada")

    original_length = int.from_bytes(data[pos:pos+8], 'big'); pos += 8
    num_ciclos = int.from_bytes(data[pos:pos+4], 'big'); pos += 4
    base = int.from_bytes(data[pos:pos+8], 'big'); pos += 8
    V = int.from_bytes(data[pos:pos+16], 'big'); pos += 16

    n_tipos = data[pos]; pos += 1
    tipos_triangulo = [data[pos+i] for i in range(n_tipos)]; pos += n_tipos

    # Cifrado
    cifrado_bytes = data[pos:]
    C = int.from_bytes(cifrado_bytes, 'big') if cifrado_bytes else 0

    # Recalcular factor
    F_total = factor_geometrico(base, tipos_triangulo, num_ciclos, V)

    # Descifrar: dividir entre el factor
    if C % F_total != 0:
        raise ValueError("El factor no divide al cifrado. Archivo corrupto o clave incorrecta.")
    M = C // F_total

    # Recuperar bytes originales (con relleno de ceros iniciales si es necesario)
    datos_originales = M.to_bytes(original_length, 'big') if original_length > 0 else b''

    with open(ruta_salida, 'wb') as f:
        f.write(datos_originales)

    print(f"✅ Archivo descifrado guardado como: {ruta_salida}")
    print(f"   Tamaño recuperado: {len(datos_originales)} bytes")

# ------------------------------------------------------------
# 4. Menú interactivo
# ------------------------------------------------------------

def menu():
    print("="*60)
    print("🔐 ENCRIPTACIÓN GEOMÉTRICA POR CONVERGENCIAS (Libro 4)")
    print("="*60)
    print("1. Encriptar archivo")
    print("2. Desencriptar archivo")
    print("3. Salir")
    return input("Elige opción (1/2/3): ").strip()

def main():
    while True:
        opcion = menu()
        if opcion == '1':
            print("\n📌 PARÁMETROS DE ENCRIPTACIÓN")
            base = int(input("  Base de la secuencia de Thales (ej: 3): ") or 3)
            tipos_input = input("  Secuencia de tipos de triángulo (ej: 1,3,2,1) [Enter para 1]: ").strip()
            if tipos_input == "":
                tipos = [1]
            else:
                tipos = [int(x) for x in tipos_input.split(',') if x.strip() != '']
            # Validar tipos
            for t in tipos:
                if t not in [1,2,3]:
                    print("❌ Tipo inválido, debe ser 1, 2 o 3.")
                    sys.exit(1)
            num_ciclos = int(input("  Número de ciclos de acumulación (ej: 10): ") or 10)

            ruta = input("\n📂 Ruta del archivo a encriptar: ").strip()
            if not os.path.isfile(ruta):
                print("❌ El archivo no existe.")
                continue
            salida = ruta + ".vma"
            encriptar_archivo(ruta, salida, base, tipos, num_ciclos)

        elif opcion == '2':
            ruta = input("\n📂 Ruta del archivo .vma a desencriptar: ").strip()
            if not os.path.isfile(ruta):
                print("❌ El archivo no existe.")
                continue
            # Quitar extensión .vma si la tiene
            if ruta.endswith('.vma'):
                salida = ruta[:-4]
            else:
                salida = ruta + ".descifrado"
            try:
                desencriptar_archivo(ruta, salida)
            except Exception as e:
                print(f"❌ Error al desencriptar: {e}")

        elif opcion == '3':
            print("👋 ¡Hasta luego!")
            break
        else:
            print("Opción no válida.")

if __name__ == "__main__":
    main()