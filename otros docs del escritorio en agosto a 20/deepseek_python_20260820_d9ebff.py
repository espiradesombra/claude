import os
import sys
import math
import hashlib
from fractions import Fraction

# ------------------------------------------------------------
# 1. Secuencias y proporciones
# ------------------------------------------------------------

def secuencia_thales(base):
    """Genera secuencia de Thales: base, base+2, base+3, base+5, base+8"""
    return [base, base+2, base+3, base+5, base+8]

RATIOS = {
    1: [Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)],          # Equilátero
    2: [Fraction(1, 4), Fraction(1, 4), Fraction(1, 2)],          # Isósceles
    3: [Fraction(1, 4), Fraction(7, 20), Fraction(2, 5)]          # Escaleno
}

# ------------------------------------------------------------
# 2. Aproximación de π con mezcla ponderada de polígonos
# ------------------------------------------------------------

def polygon_pi(sides, iterations):
    """
    Aproxima π usando un polígono regular de 'sides' lados.
    Se parte de un ángulo inicial y se duplica el número de lados
    'iterations' veces, acumulando el perímetro.
    Devuelve la estimación de π.
    """
    # Ángulo central inicial (usa math.pi, pero solo como semilla educativa)
    angulo = math.pi / sides
    lado = 2 * math.sin(angulo)
    n = sides
    for _ in range(iterations):
        # Fórmula de duplicación del lado en polígono inscrito
        lado = lado / math.sqrt(2 + math.sqrt(4 - lado**2))
        n *= 2
    perimetro = n * lado
    return perimetro / 2

def weighted_pi(components):
    """
    components: lista de tuplas (peso%, lados, iteraciones)
    Ej: [(30,3,70), (20,7,66), (50,4,444)]
    Devuelve la aproximación ponderada de π.
    """
    total_peso = sum(w for w, _, _ in components)
    suma = 0.0
    for peso, lados, iters in components:
        pi_i = polygon_pi(lados, iters)
        suma += peso * pi_i
    return suma / total_peso

def approx_e(iterations=50):
    """Aproxima e usando serie de Taylor con 'iterations' términos."""
    e = 0.0
    factorial = 1.0
    for k in range(iterations):
        if k == 0:
            e += 1.0
        else:
            factorial *= k
            e += 1.0 / factorial
    return e

# ------------------------------------------------------------
# 3. Derivación determinista del factor geométrico
# ------------------------------------------------------------

def derive_v(base, tipos, ciclos, pi_aprox, e_aprox):
    """
    Deriva un entero V a partir de todos los parámetros (incluidos π y e).
    Así no se guarda V en los metadatos; se recalcula al desencriptar.
    """
    # Formateamos π y e con precisión fija para evitar ambigüedades
    pi_str = f"{pi_aprox:.15e}"
    e_str = f"{e_aprox:.15e}"
    cadena = f"{base}|{tipos}|{ciclos}|{pi_str}|{e_str}"
    digest = hashlib.sha256(cadena.encode('utf-8')).digest()
    V = int.from_bytes(digest[:16], 'big') + 1
    return V

def factor_geometrico(base, tipos, ciclos, pi_aprox, e_aprox):
    """
    Calcula el factor (1 + a) * V, donde a es la acumulación de lados
    según la secuencia de Thales y los tipos de triángulo.
    Devuelve el factor entero F_total.
    """
    S = secuencia_thales(base)
    a = Fraction(0)
    for i in range(ciclos):
        tipo = tipos[i % len(tipos)]
        escala = S[i % len(S)]
        for r in RATIOS[tipo]:
            a += r * escala
    a = a % 1  # parte fraccionaria
    factor_frac = 1 + a
    F = factor_frac.numerator

    V = derive_v(base, tipos, ciclos, pi_aprox, e_aprox)
    return F * V

# ------------------------------------------------------------
# 4. Encriptar / desencriptar
# ------------------------------------------------------------

def encriptar_archivo(ruta_entrada, ruta_salida, base, tipos, ciclos,
                      pi_components, e_iters, guardar_metadatos=True):
    # Leer archivo
    with open(ruta_entrada, 'rb') as f:
        datos = f.read()

    original_length = len(datos)
    M = int.from_bytes(datos, 'big') if datos else 0

    # Calcular π y e
    pi_aprox = weighted_pi(pi_components)
    e_aprox = approx_e(e_iters)

    # Factor geométrico
    F_total = factor_geometrico(base, tipos, ciclos, pi_aprox, e_aprox)

    # Cifrar: multiplicar el segmento por el factor
    C = M * F_total

    # Convertir a bytes
    cifrado_bytes = C.to_bytes((C.bit_length() + 7) // 8, 'big') if C > 0 else b'\x00'

    # Construir encabezado
    magic = b'GEO2'
    version = 1
    flags = 0
    if guardar_metadatos:
        flags |= 0b00000001  # bit0: metadatos presentes

    header = magic + bytes([version, flags])
    header += original_length.to_bytes(8, 'big')

    if guardar_metadatos:
        header += base.to_bytes(8, 'big')
        header += bytes([len(tipos)])
        for t in tipos:
            header += bytes([t])
        header += ciclos.to_bytes(4, 'big')
        # Guardamos π y e como cadenas de longitud fija
        pi_bytes = f"{pi_aprox:.15e}".encode('utf-8')
        e_bytes = f"{e_aprox:.15e}".encode('utf-8')
        header += bytes([len(pi_bytes)]) + pi_bytes
        header += bytes([len(e_bytes)]) + e_bytes

    with open(ruta_salida, 'wb') as f:
        f.write(header + cifrado_bytes)

    print(f"✅ Archivo encriptado guardado como: {ruta_salida}")
    print(f"   Tamaño original: {original_length} bytes")
    print(f"   Tamaño cifrado total: {len(header) + len(cifrado_bytes)} bytes")
    print(f"   π aprox: {pi_aprox:.15f}")
    print(f"   e aprox: {e_aprox:.15f}")
    print(f"   Factor F_total: {F_total}")

def desencriptar_core(cifrado_bytes, original_length, base, tipos, ciclos,
                      pi_aprox, e_aprox):
    C = int.from_bytes(cifrado_bytes, 'big') if cifrado_bytes else 0
    F_total = factor_geometrico(base, tipos, ciclos, pi_aprox, e_aprox)

    if C % F_total != 0:
        raise ValueError("El factor no divide al cifrado. Parámetros incorrectos.")

    M = C // F_total
    datos = M.to_bytes(original_length, 'big') if original_length > 0 else b''
    return datos

def desencriptar_archivo(ruta_entrada, ruta_salida, params_manual=None):
    with open(ruta_entrada, 'rb') as f:
        data = f.read()

    # Parsear cabecera
    magic = data[0:4]
    if magic != b'GEO2':
        raise ValueError("Formato no reconocido")

    version = data[4]
    flags = data[5]
    pos = 6
    original_length = int.from_bytes(data[pos:pos+8], 'big'); pos += 8

    if flags & 0b00000001:
        # Metadatos presentes
        base = int.from_bytes(data[pos:pos+8], 'big'); pos += 8
        n_tipos = data[pos]; pos += 1
        tipos = [data[pos+i] for i in range(n_tipos)]; pos += n_tipos
        ciclos = int.from_bytes(data[pos:pos+4], 'big'); pos += 4

        len_pi = data[pos]; pos += 1
        pi_bytes = data[pos:pos+len_pi]; pos += len_pi
        len_e = data[pos]; pos += 1
        e_bytes = data[pos:pos+len_e]; pos += len_e
        pi_aprox = float(pi_bytes.decode('utf-8'))
        e_aprox = float(e_bytes.decode('utf-8'))
    else:
        # Sin metadatos: el usuario debe pasar los parámetros manualmente
        if params_manual is None:
            raise ValueError("Este archivo no tiene metadatos. Proporcione los parámetros.")
        base = params_manual['base']
        tipos = params_manual['tipos']
        ciclos = params_manual['ciclos']
        pi_components = params_manual['pi_components']
        e_iters = params_manual['e_iters']
        pi_aprox = weighted_pi(pi_components)
        e_aprox = approx_e(e_iters)

    cifrado_bytes = data[pos:]

    datos = desencriptar_core(cifrado_bytes, original_length, base, tipos,
                              ciclos, pi_aprox, e_aprox)

    with open(ruta_salida, 'wb') as f:
        f.write(datos)

    print(f"✅ Archivo descifrado guardado como: {ruta_salida}")
    print(f"   Tamaño recuperado: {len(datos)} bytes")

# ------------------------------------------------------------
# 5. Menú interactivo educativo
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
            print("\n📌 PARÁMETROS DE ENCRIPTACIÓN (modo educativo)")
            print("   La base genera una secuencia de Thales para escalar los triángulos.")
            base = int(input("   Base de la secuencia de Thales (ej: 3): ") or 3)

            print("   Introduce una secuencia de tipos de triángulo (1=Equilátero, 2=Isósceles, 3=Escaleno).")
            print("   Se usará cíclicamente para acumular lados.")
            tipos_input = input("   Secuencia (ej: 1,3,2,1) [Enter para 1]: ").strip()
            if tipos_input == "":
                tipos = [1]
            else:
                tipos = [int(x) for x in tipos_input.split(',') if x.strip() != '']
            for t in tipos:
                if t not in [1,2,3]:
                    print("❌ Tipo inválido, debe ser 1, 2 o 3.")
                    sys.exit(1)

            print("   Número de ciclos: cuántas veces se recorre la secuencia de triángulos.")
            ciclos = int(input("   Ciclos de acumulación (ej: 3): ") or 3)

            print("\n🧮 APROXIMACIÓN DE π POR MEZCLA DE POLÍGONOS")
            print("   Puedes combinar varios polígonos con pesos (%).")
            print("   Ejemplo: 30% triángulo 70 iter, 20% heptágono 66 iter, 50% cuadrado 444 iter.")
            n_comp = int(input("   Número de componentes para π (ej: 3): ") or 1)
            pi_components = []
            for i in range(n_comp):
                print(f"   Componente {i+1}:")
                peso = float(input("     Peso % (ej: 30): ") or 30)
                lados = int(input("     Nº lados del polígono (ej: 3): ") or 3)
                iters = int(input("     Iteraciones de duplicación (ej: 70): ") or 70)
                pi_components.append((peso, lados, iters))

            print("\n📐 APROXIMACIÓN DE e")
            e_iters = int(input("   Iteraciones de la serie de Taylor para e (ej: 50): ") or 50)

            print("\n💾 METADATOS")
            guardar = input("   ¿Guardar parámetros dentro del archivo .vma? (s/n) [s]: ").strip().lower()
            guardar_metadatos = (guardar != 'n')

            ruta = input("\n📂 Ruta del archivo a encriptar: ").strip()
            if not os.path.isfile(ruta):
                print("❌ El archivo no existe.")
                continue
            salida = ruta + ".vma"
            encriptar_archivo(ruta, salida, base, tipos, ciclos,
                              pi_components, e_iters, guardar_metadatos)

        elif opcion == '2':
            ruta = input("\n📂 Ruta del archivo .vma a desencriptar: ").strip()
            if not os.path.isfile(ruta):
                print("❌ El archivo no existe.")
                continue
            if ruta.endswith('.vma'):
                salida = ruta[:-4]
            else:
                salida = ruta + ".descifrado"

            try:
                # Intentamos desencriptar; si no hay metadatos, se piden manualmente
                desencriptar_archivo(ruta, salida)
            except ValueError as e:
                if "sin metadatos" in str(e):
                    print("\n⚠️  El archivo no incluye metadatos. Debes introducir los mismos parámetros que usaste al encriptar.")
                    base = int(input("   Base de secuencia de Thales: "))
                    tipos_input = input("   Secuencia de tipos de triángulo (ej: 1,3,2,1): ").strip()
                    tipos = [int(x) for x in tipos_input.split(',') if x.strip() != '']
                    ciclos = int(input("   Ciclos de acumulación: "))
                    n_comp = int(input("   Número de componentes para π: "))
                    pi_components = []
                    for i in range(n_comp):
                        peso = float(input(f"     Componente {i+1} peso %: "))
                        lados = int(input(f"     Componente {i+1} lados: "))
                        iters = int(input(f"     Componente {i+1} iteraciones: "))
                        pi_components.append((peso, lados, iters))
                    e_iters = int(input("   Iteraciones para e: "))
                    params_manual = {
                        'base': base,
                        'tipos': tipos,
                        'ciclos': ciclos,
                        'pi_components': pi_components,
                        'e_iters': e_iters
                    }
                    desencriptar_archivo(ruta, salida, params_manual)
                else:
                    print(f"❌ Error: {e}")
        elif opcion == '3':
            print("👋 ¡Hasta luego!")
            break
        else:
            print("Opción no válida.")

if __name__ == "__main__":
    main()