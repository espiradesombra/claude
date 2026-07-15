"""Criba Desmemoriada — primos 6k±1 (Python, reproducible en Windows)."""


def criba_desmemoriada(limite: int) -> list[int]:
    p6 = [True] * (limite + 1)
    p6[0] = p6[1] = False
    for i in range(2, limite + 1):
        if i % 2 == 0 or i % 3 == 0:
            p6[i] = False
    for p in range(5, int(limite**0.5) + 1, 6):
        if p6[p]:
            for multiple in range(p * p, limite + 1, 6 * p):
                p6[multiple] = False
                if multiple + 2 * p <= limite:
                    p6[multiple + 2 * p] = False
    return [i for i, is_prime in enumerate(p6) if is_prime]


if __name__ == "__main__":
    primos = criba_desmemoriada(100)
    print(f"Primos hasta 100: {len(primos)}")
    print(primos)