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
    """
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

# ------------------------------------------------------------
# 3. Aproximación de e personalizable
# ------------------------------------------------------------

def approx_e_taylor(iterations=50):
    """Aproxima e usando serie de Taylor clásica."""
    e = 0.0
    factorial = 1.0
    for k in range(iterations):
        if k == 0:
            e += 1.0
        else:
            factorial *= k
            e += 1.0 / factorial
    return e

def approx_e_custom(skip_modulos, n_terms):
    """
    Aproximación de e sumando 1/n, descartando múltiplos de skip_modulos.
    Al final se suma 2.
    """
    acum = 0.0
    count = 0
    i = 1
    while count < n_terms:
        # Si i es múltiplo de alguno de los módulos prohibidos, se omite
        if any(i % m == 0 for m in skip_modulos):
            i += 1
            continue
        acum += 1.0 / i
        count += 1
        i += 1
    return 2.0 + acum

def compute_e_from_params(e_params):
    """
    e_params: dict con modo y parámetros.
    Modos:
      - 'taylor': {'mode':'taylor', 'iters': N}
      - 'custom': {'mode':'custom', 'mods': lista, 'terms': N}
      - 'double': {'mode':'double', 'mods1': lista, 'terms1': N,
                   'mods2': lista, 'terms2': M}
    """
    mode = e_params.get('mode', 'taylor')
    if mode == 'taylor':
        return approx_e_taylor(e_params.get('iters', 50))
    elif mode == 'custom':
        return approx_e_custom(e_params.get('mods', [2]), e_params.get('terms', 50))
    elif mode == 'double':
        e1 = approx_e_custom(e_params.get('mods1', [2]), e_params.get('terms1', 50))
        e2 = approx_e_custom(e_params.get('mods2', [7]), e_params.get('terms2', 60))
        return (e1 + e2) / 2.0
    else:
        return approx_e_taylor(50)

# ------------------------------------------------------------
# 4. Derivación determinista del factor geométrico
# ------------------------------------------------------------

def derive_v(base, tipos, ciclos, pi_str, e_str, omitir_patron):
    """
    Deriva un entero V a partir de todos los parámetros.
    pi_str y e_str son cadenas con formato estable.
    """
    omitir_str = ''.join(str(x) for x in omitir_patron)
    cadena = f"{base}|{tipos}|{ciclos}|{pi_str}|{e_str}|{omitir_str}"
    digest = hashlib.sha256(cadena.encode('utf-8')).digest()
    V = int.from_bytes(digest[:16], 'big') + 1
    return V

def factor_geometrico(base, tipos, ciclos, pi_str, e_str, omitir_patron):
    """
    Calcula el factor (1 + a) * V, donde a es la acumulación de lados
    según la secuencia de Thales, los tipos de triángulo y los lados omitidos.
    """
    S = secuencia_thales(base)
    a = Fraction(0)

    for i in range(ciclos):
        tipo = tipos[i % len(tipos)]
        escala = S[i % len(S)]
        ratios = RATIOS[tipo]
        for j, r in enumerate(ratios):
            # Índice global del lado (cada triángulo tiene 3 lados)
            idx_global = i * 3 + j
            # Si el patrón indica 0, se omite este lado
            if omitir_patron[idx_global % len(omitir_patron)] == 0:
                continue
            a += r * escala

    a = a % 1                     # parte fraccionaria
    factor_frac = 1 + a
    F = factor_frac.numerator

    V = derive_v(base, tipos, ciclos, pi_str, e_str, omitir_patron)
    return F * V

# ------------------------------------------------------------
# 5. Encriptar / desencriptar
# ------------------------------------------------------------

def encriptar_archivo(ruta_entrada, ruta_salida, base, tipos, ciclos,
                      pi_components, e_params, omitir_patron, guardar_metadatos=True):
    # Leer archivo
    with open(ruta_entrada, 'rb') as f:
        datos = f.read()

    original_length = len(datos)
    M = int.from_bytes(datos, 'big') if datos else 0

    # Calcular π y e como flotantes y luego como cadenas estables
    pi_aprox = weighted_pi(pi_components)
    e_aprox = compute_e_from_params(e_params)
    pi_str = f"{pi_aprox:.15e}"
    e_str = f"{e_aprox:.15e}"

    # Factor geométrico
    F_total = factor_geometrico(base, tipos, ciclos, pi_str, e_str, omitir_patron)

    # Cifrar: multiplicar el segmento por el factor
    C = M * F_total

    # Convertir a bytes
    cifrado_bytes = C.to_bytes((C.bit_length() + 7) // 8, 'big') if C > 0 else b'\x00'

    # Construir encabezado
    magic = b'GEO3'
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
        # Guardar patrón de omisión
        header += bytes([len(omitir_patron)])
        for v in omitir_patron:
            header += bytes([v])
        # Guardar π y e como cadenas de longitud fija
        pi_bytes = pi_str.encode('utf-8')
        e_bytes = e_str.encode('utf-8')
        header += bytes([len(pi_bytes)]) + pi_bytes
        header += bytes([len(e_bytes)]) + e_bytes

    with open(ruta_salida, 'wb') as f:
        f.write(header + cifrado_bytes)

    print(f"✅ Archivo encriptado guardado como: {ruta_salida}")
    print(f"   Tamaño original: {original_length} bytes")
    print(f"   Tamaño cifrado total: {len(header) + len(cifrado_bytes)} bytes")
    print(f"   π aprox: {pi_aprox:.15f}")
    print(f"   e aprox: {e_aprox:.15f}")
    print(f"   Patrón de lados omitidos: {omitir_patron}")
    print(f"   Factor F_total: {F_total}")

def desencriptar_core(cifrado_bytes, original_length, base, tipos, ciclos,
                      pi_str, e_str, omitir_patron):
    C = int.from_bytes(cifrado_bytes, 'big') if cifrado_bytes else 0
    F_total = factor_geometrico(base, tipos, ciclos, pi_str, e_str, omitir_patron)

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
    if magic != b'GEO3':
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

        n_omitir = data[pos]; pos += 1
        omitir_patron = [data[pos+i] for i in range(n_omitir)]; pos += n_omitir

        len_pi = data[pos]; pos += 1
        pi_bytes = data[pos:pos+len_pi]; pos += len_pi
        len_e = data[pos]; pos += 1
        e_bytes = data[pos:pos+len_e]; pos += len_e
        pi_str = pi_bytes.decode('utf-8')
        e_str = e_bytes.decode('utf-8')
    else:
        # Sin metadatos: el usuario debe pasar los parámetros manualmente
        if params_manual is None:
            raise ValueError("Este archivo no tiene metadatos. Proporcione los parámetros.")
        base = params_manual['base']
        tipos = params_manual['tipos']
        ciclos = params_manual['ciclos']
        omitir_patron = params_manual['omitir_patron']
        pi_components = params_manual['pi_components']
        e_params = params_manual['e_params']
        pi_aprox = weighted_pi(pi_components)
        e_aprox = compute_e_from_params(e_params)
        pi_str = f"{pi_aprox:.15e}"
        e_str = f"{e_aprox:.15e}"

    cifrado_bytes = data[pos:]

    datos = desencriptar_core(cifrado_bytes, original_length, base, tipos,
                              ciclos, pi_str, e_str, omitir_patron)

    with open(ruta_salida, 'wb') as f:
        f.write(datos)

    print(f"✅ Archivo descifrado guardado como: {ruta_salida}")
    print(f"   Tamaño recuperado: {len(datos)} bytes")

# ------------------------------------------------------------
# 6. Menú interactivo educativo
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

            print("\n   Introduce una secuencia de tipos de triángulo (1=Equilátero, 2=Isósceles, 3=Escaleno).")
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

            print("\n   Número de ciclos: cuántas veces se recorre la secuencia de triángulos.")
            ciclos = int(input("   Ciclos de acumulación (ej: 3): ") or 3)

            print("\n🔘 Secuencia de lados que NO cuentan")
            print("   Introduce un patrón de 0 y 1 separados por comas.")
            print("   1 = lado cuenta, 0 = lado omitido.")
            omitir_input = input("   Patrón (ej: 1,0,1) [Enter para 1,1,1]: ").strip()
            if omitir_input == "":
                omitir_patron = [1,1,1]
            else:
                omitir_patron = [int(x) for x in omitir_input.split(',') if x.strip() != '']
            # Validar que sean 0 o 1
            for v in omitir_patron:
                if v not in [0,1]:
                    print("❌ El patrón debe contener solo 0 o 1.")
                    sys.exit(1)

            print("\n🧮 APROXIMACIÓN DE π POR MEZCLA DE POLÍGONOS")
            print("   Puedes combinar varios polígonos con pesos (%).")
            n_comp = int(input("   Número de componentes para π (ej: 2): ") or 1)
            pi_components = []
            for i in range(n_comp):
                print(f"   Componente {i+1}:")
                peso = float(input("     Peso % (ej: 40): ") or 40)
                lados = int(input("     Nº lados del polígono (ej: 6): ") or 6)
                iters = int(input("     Iteraciones de duplicación (ej: 66): ") or 66)
                pi_components.append((peso, lados, iters))

            print("\n📐 APROXIMACIÓN DE e")
            print("   Elige método:")
            print("   1: Serie de Taylor clásica")
            print("   2: Descarte de múltiplos (custom)")
            print("   3: Descarte doble y media")
            modo_e = int(input("   Método (1/2/3): ") or 1)
            if modo_e == 1:
                e_iters = int(input("   Iteraciones Taylor (ej: 50): ") or 50)
                e_params = {'mode': 'taylor', 'iters': e_iters}
            elif modo_e == 2:
                mods = input("   Lista de módulos a descartar (ej: 2,3,5): ").strip()
                mods = [int(x) for x in mods.split(',') if x.strip()] if mods else [2]
                terms = int(input("   Número de términos válidos (ej: 60): ") or 60)
                e_params = {'mode': 'custom', 'mods': mods, 'terms': terms}
            elif modo_e == 3:
                mods1 = input("   Módulos primera pasada (ej: 2,3,5): ").strip()
                mods1 = [int(x) for x in mods1.split(',') if x.strip()] if mods1 else [2]
                terms1 = int(input("   Términos primera pasada (ej: 50): ") or 50)
                mods2 = input("   Módulos segunda pasada (ej: 7,11): ").strip()
                mods2 = [int(x) for x in mods2.split(',') if x.strip()] if mods2 else [7]
                terms2 = int(input("   Términos segunda pasada (ej: 60): ") or 60)
                e_params = {'mode': 'double', 'mods1': mods1, 'terms1': terms1,
                            'mods2': mods2, 'terms2': terms2}
            else:
                print("❌ Método inválido. Se usará Taylor.")
                e_params = {'mode': 'taylor', 'iters': 50}

            print("\n💾 METADATOS")
            guardar = input("   ¿Guardar parámetros dentro del archivo .vma? (s/n) [s]: ").strip().lower()
            guardar_metadatos = (guardar != 'n')

            ruta = input("\n📂 Ruta del archivo a encriptar: ").strip()
            if not os.path.isfile(ruta):
                print("❌ El archivo no existe.")
                continue
            salida = ruta + ".vma"
            encriptar_archivo(ruta, salida, base, tipos, ciclos,
                              pi_components, e_params, omitir_patron, guardar_metadatos)

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
                desencriptar_archivo(ruta, salida)
            except ValueError as e:
                if "sin metadatos" in str(e):
                    print("\n⚠️  El archivo no incluye metadatos. Debes introducir los mismos parámetros que usaste al encriptar.")
                    base = int(input("   Base de secuencia de Thales: "))
                    tipos_input = input("   Secuencia de tipos de triángulo (ej: 1,3,2,1): ").strip()
                    tipos = [int(x) for x in tipos_input.split(',') if x.strip() != '']
                    ciclos = int(input("   Ciclos de acumulación: "))
                    omitir_input = input("   Patrón de lados omitidos (ej: 1,0,1): ").strip()
                    omitir_patron = [int(x) for x in omitir_input.split(',') if x.strip() != '']
                    n_comp = int(input("   Número de componentes para π: "))
                    pi_components = []
                    for i in range(n_comp):
                        peso = float(input(f"     Componente {i+1} peso %: "))
                        lados = int(input(f"     Componente {i+1} lados: "))
                        iters = int(input(f"     Componente {i+1} iteraciones: "))
                        pi_components.append((peso, lados, iters))
                    print("   Parámetros para e:")
                    modo_e = int(input("   Método e (1=Taylor, 2=Custom, 3=Double): ") or 1)
                    if modo_e == 1:
                        e_iters = int(input("   Iteraciones Taylor: ") or 50)
                        e_params = {'mode': 'taylor', 'iters': e_iters}
                    elif modo_e == 2:
                        mods = input("   Módulos a descartar: ").strip()
                        mods = [int(x) for x in mods.split(',') if x.strip()] if mods else [2]
                        terms = int(input("   Términos: ") or 60)
                        e_params = {'mode': 'custom', 'mods': mods, 'terms': terms}
                    else:
                        mods1 = input("   Módulos primera pasada: ").strip()
                        mods1 = [int(x) for x in mods1.split(',') if x.strip()] if mods1 else [2]
                        terms1 = int(input("   Términos primera pasada: ") or 50)
                        mods2 = input("   Módulos segunda pasada: ").strip()
                        mods2 = [int(x) for x in mods2.split(',') if x.strip()] if mods2 else [7]
                        terms2 = int(input("   Términos segunda pasada: ") or 60)
                        e_params = {'mode': 'double', 'mods1': mods1, 'terms1': terms1,
                                    'mods2': mods2, 'terms2': terms2}

                    params_manual = {
                        'base': base,
                        'tipos': tipos,
                        'ciclos': ciclos,
                        'omitir_patron': omitir_patron,
                        'pi_components': pi_components,
                        'e_params': e_params
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