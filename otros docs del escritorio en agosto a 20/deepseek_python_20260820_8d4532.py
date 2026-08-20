import os
import sys
import math
import hashlib
import hmac
import secrets
import getpass

# ------------------------------------------------------------
# 1. Aproximaciones de π y e (según el Libro 4)
# ------------------------------------------------------------

def aprox_pi(n_lados=100000, iteraciones=5):
    """
    Aproxima π usando un polígono regular de n_lados con duplicación.
    """
    # Usamos math.pi como valor inicial para el ángulo (simplificado)
    angulo = math.pi / n_lados
    lado = 2 * math.sin(angulo)
    perimetro = n_lados * lado
    pi_aprox = perimetro / 2

    for _ in range(iteraciones):
        # Fórmula de duplicación del lado
        lado = lado / math.sqrt(2 + math.sqrt(4 - lado**2))
        n_lados *= 2
        perimetro = n_lados * lado
        pi_aprox = perimetro / 2

    return pi_aprox

def aprox_e(iteraciones=50):
    """
    Aproxima e usando serie de Taylor.
    """
    e_aprox = 0.0
    factorial = 1.0
    for k in range(iteraciones):
        if k == 0:
            e_aprox += 1.0
        else:
            factorial *= k
            e_aprox += 1.0 / factorial
    return e_aprox

# ------------------------------------------------------------
# 2. Secuencia de Thales y tipos de triángulo
# ------------------------------------------------------------

def secuencia_thales(base):
    """Genera la secuencia de Thales a partir de una base."""
    return [base, base+2, base+3, base+5, base+8]

def tipo_figura(opcion):
    """Devuelve las proporciones de lados del triángulo seleccionado."""
    if opcion == 1:
        return [1/3, 1/3, 1/3]      # Equilátero
    elif opcion == 2:
        return [1/4, 1/4, 1/2]      # Isósceles
    elif opcion == 3:
        return [0.25, 0.35, 0.40]   # Escaleno
    else:
        return [1/3, 1/3, 1/3]

# ------------------------------------------------------------
# 3. Derivación de clave a partir de los parámetros
# ------------------------------------------------------------

def derivar_clave(params, salt):
    """
    Deriva una clave de 32 bytes usando SHA-256 sobre los parámetros y la sal.
    """
    # Convertir todos los parámetros a una cadena canónica
    cadena = f"{params['secuencia']}|{params['figura']}|{params['iteraciones']}|" \
             f"{params['pi_aprox']:.15f}|{params['e_aprox']:.15f}|{salt.hex()}"
    # Hash SHA-256
    return hashlib.sha256(cadena.encode('utf-8')).digest()

# ------------------------------------------------------------
# 4. Cifrado / descifrado con keystream y HMAC
# ------------------------------------------------------------

def cifrar_bytes(datos, clave, salt):
    """
    Cifra usando XOR con un keystream generado por SHA-256 en contador.
    """
    # Inicializamos el contador
    contador = 0
    keystream = bytearray()
    while len(keystream) < len(datos):
        # Bloque de keystream: SHA-256(salt || contador)
        bloque_entrada = salt + contador.to_bytes(8, 'big')
        bloque_salida = hashlib.sha256(bloque_entrada).digest()
        keystream.extend(bloque_salida)
        contador += 1
    # XOR
    datos_cifrados = bytes(b ^ k for b, k in zip(datos, keystream))
    # HMAC para integridad
    mac = hmac.new(clave, datos_cifrados, hashlib.sha256).digest()
    return datos_cifrados, mac

def descifrar_bytes(datos_cifrados, clave, salt, mac_recibido):
    """
    Descifra verificando HMAC.
    """
    # Verificar HMAC
    mac_calculado = hmac.new(clave, datos_cifrados, hashlib.sha256).digest()
    if not hmac.compare_digest(mac_calculado, mac_recibido):
        raise ValueError("❌ Error de integridad: el archivo fue alterado o la clave es incorrecta.")
    # Regenerar keystream igual que en cifrado
    contador = 0
    keystream = bytearray()
    while len(keystream) < len(datos_cifrados):
        bloque_entrada = salt + contador.to_bytes(8, 'big')
        bloque_salida = hashlib.sha256(bloque_entrada).digest()
        keystream.extend(bloque_salida)
        contador += 1
    datos = bytes(b ^ k for b, k in zip(datos_cifrados, keystream))
    return datos

# ------------------------------------------------------------
# 5. Funciones de archivo
# ------------------------------------------------------------

def encriptar_archivo(ruta_entrada, ruta_salida, params):
    """
    Lee el archivo, lo cifra y guarda con formato .vma
    """
    # Leer archivo original
    with open(ruta_entrada, 'rb') as f:
        datos = f.read()

    # Generar sal aleatoria (16 bytes)
    salt = secrets.token_bytes(16)

    # Derivar clave
    clave = derivar_clave(params, salt)

    # Cifrar
    datos_cifrados, mac = cifrar_bytes(datos, clave, salt)

    # Guardar archivo cifrado con metadatos
    # Formato: [magic(4)] [salt(16)] [mac(32)] [longitud_params(2)] [params_json] [datos_cifrados]
    magic = b'VMA1'
    import json
    params_json = json.dumps({
        'secuencia': params['secuencia'],
        'figura': params['figura'],
        'iteraciones': params['iteraciones'],
        'pi_aprox': params['pi_aprox'],
        'e_aprox': params['e_aprox']
    }).encode('utf-8')
    longitud_params = len(params_json).to_bytes(2, 'big')

    with open(ruta_salida, 'wb') as f:
        f.write(magic)
        f.write(salt)
        f.write(mac)
        f.write(longitud_params)
        f.write(params_json)
        f.write(datos_cifrados)

    print(f"✅ Archivo encriptado guardado como: {ruta_salida}")
    print(f"   Tamaño original: {len(datos)} bytes")
    print(f"   Tamaño cifrado: {len(datos_cifrados)} bytes")

def desencriptar_archivo(ruta_entrada, ruta_salida):
    """
    Lee un archivo .vma, verifica integridad y lo descifra.
    """
    with open(ruta_entrada, 'rb') as f:
        magic = f.read(4)
        if magic != b'VMA1':
            raise ValueError("❌ Formato no reconocido (falta magic 'VMA1')")
        salt = f.read(16)
        mac_recibido = f.read(32)
        longitud_params = int.from_bytes(f.read(2), 'big')
        params_json = f.read(longitud_params)
        datos_cifrados = f.read()

    # Parsear parámetros
    import json
    params_dict = json.loads(params_json.decode('utf-8'))
    params = {
        'secuencia': params_dict['secuencia'],
        'figura': params_dict['figura'],
        'iteraciones': params_dict['iteraciones'],
        'pi_aprox': params_dict['pi_aprox'],
        'e_aprox': params_dict['e_aprox']
    }

    # Derivar clave
    clave = derivar_clave(params, salt)

    # Descifrar (verifica HMAC)
    datos = descifrar_bytes(datos_cifrados, clave, salt, mac_recibido)

    # Guardar archivo descifrado
    with open(ruta_salida, 'wb') as f:
        f.write(datos)

    print(f"✅ Archivo descifrado guardado como: {ruta_salida}")
    print(f"   Tamaño recuperado: {len(datos)} bytes")

# ------------------------------------------------------------
# 6. Interfaz de menú para Colab
# ------------------------------------------------------------

def menu():
    print("="*60)
    print("🔐 ENCRIPTACIÓN GEOMÉTRICA — VERSIÓN COLAB (Libro 4)")
    print("="*60)
    print("1. Encriptar archivo")
    print("2. Desencriptar archivo")
    print("3. Salir")
    opcion = input("Elige opción (1/2/3): ").strip()
    return opcion

def main():
    while True:
        opcion = menu()
        if opcion == '1':
            # --- ENCRIPTAR ---
            print("\n📌 PARÁMETROS DE ENCRIPTACIÓN")
            base = int(input("  Base de secuencia de Thales (ej: 3): ") or 3)
            secuencia = secuencia_thales(base)
            print(f"    Secuencia: {secuencia}")

            print("  Tipos de figura: 1=Equilátero, 2=Isósceles, 3=Escaleno")
            tipo = int(input("  Elige tipo (1/2/3): ") or 1)
            figura = tipo_figura(tipo)
            print(f"    Figura: {figura}")

            iteraciones = int(input("  Iteraciones de ofuscación (ej: 10): ") or 10)

            n_lados = int(input("  Nº lados para π (ej: 100000): ") or 100000)
            iter_pi = int(input("  Iteraciones duplicación π (ej: 5): ") or 5)
            iter_e = int(input("  Iteraciones serie e (ej: 50): ") or 50)

            print("\n⏳ Calculando aproximaciones...")
            pi_aprox = aprox_pi(n_lados, iter_pi)
            e_aprox = aprox_e(iter_e)
            print(f"    π ≈ {pi_aprox:.15f}")
            print(f"    e ≈ {e_aprox:.15f}")

            params = {
                'secuencia': secuencia,
                'figura': figura,
                'iteraciones': iteraciones,
                'pi_aprox': pi_aprox,
                'e_aprox': e_aprox
            }

            # Ruta del archivo (en Colab)
            ruta = input("\n📂 Ruta del archivo a encriptar: ").strip()
            if not os.path.isfile(ruta):
                print("❌ El archivo no existe.")
                continue
            salida = ruta + ".vma"
            encriptar_archivo(ruta, salida, params)

        elif opcion == '2':
            # --- DESENCRIPTAR ---
            ruta = input("\n📂 Ruta del archivo .vma a desencriptar: ").strip()
            if not os.path.isfile(ruta):
                print("❌ El archivo no existe.")
                continue
            # Quitar extensión .vma si la tiene
            salida = ruta
            if salida.endswith('.vma'):
                salida = salida[:-4]
            else:
                salida = salida + ".descifrado"
            try:
                desencriptar_archivo(ruta, salida)
            except Exception as e:
                print(f"❌ Error al desencriptar: {e}")

        elif opcion == '3':
            print("👋 ¡Hasta luego!")
            break
        else:
            print("Opción no válida.")

# Ejecutar menú
if __name__ == "__main__":
    main()