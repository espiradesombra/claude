# Dentro de main(), en la opción de encriptar:
print("\n📌 PARÁMETROS DE ENCRIPTACIÓN")
base = int(input("  Base de secuencia de Thales (ej: 3): ") or 3)
secuencia = secuencia_thales(base)
print(f"    Secuencia Thales: {secuencia}")

# --- NUEVO: secuencia de tipos de triángulo ---
tipos_input = input("  Secuencia de tipos de triángulo (ej: 1,3,2,1) [Enter para 1]: ").strip()
if tipos_input == "":
    tipos_triangulo = [1]  # equilátero fijo
else:
    tipos_triangulo = [int(x) for x in tipos_input.split(',') if x.strip() != '']
# Validar que todos sean 1, 2 o 3
for t in tipos_triangulo:
    if t not in [1,2,3]:
        print("❌ Tipo inválido, debe ser 1, 2 o 3.")
        sys.exit(1)
print(f"    Secuencia de tipos: {tipos_triangulo}")

iteraciones = int(input("  Iteraciones de ofuscación (ej: 10): ") or 10)
n_lados = int(input("  Nº lados para π (ej: 100000): ") or 100000)
iter_pi = int(input("  Iteraciones duplicación π (ej: 5): ") or 5)
iter_e = int(input("  Iteraciones serie e (ej: 50): ") or 50)