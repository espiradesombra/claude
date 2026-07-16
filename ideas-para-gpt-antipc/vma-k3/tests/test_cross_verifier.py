import unittest

from vma.k3 import K3CrossVerifier, K3_FIRMA_LOGICA
from vma.k3.toffoli_hash import theorem_k3_hash


class TestK3CrossVerifier(unittest.TestCase):
    def test_firma_logica_presente(self) -> None:
        self.assertIn(b"9E3779B1", K3_FIRMA_LOGICA)

    def test_deterministic(self) -> None:
        data = b"LECTURA_DE_SENSORES_DATA"
        v = K3CrossVerifier(bits_ancho=32, desfase_stride=0, num_registros=3, semilla=0x1F2E3D4C)
        h1 = v.verificar_memoria(data)
        h2 = v.verificar_memoria(data)
        self.assertEqual(h1, h2)
        self.assertRegex(h1, r"^0x[0-9A-F]{8}$")

    def test_stride_changes_hash(self) -> None:
        data = b"\xAA\xBB\x12\x34\xCC\xDD\x56\x78"
        v0 = K3CrossVerifier(bits_ancho=16, desfase_stride=0, num_registros=3, semilla=0x1F2E3D4C)
        v2 = K3CrossVerifier(bits_ancho=16, desfase_stride=2, num_registros=3, semilla=0x1F2E3D4C)
        self.assertNotEqual(v0.verificar_memoria(data), v2.verificar_memoria(data))

    def test_marcas_change_hash(self) -> None:
        data = b"test"
        v = K3CrossVerifier(semilla=0x1F2E3D4C)
        base = v.verificar_memoria(data)
        v.encolar_marca_final(0x04, 0x10)
        self.assertNotEqual(base, v.verificar_memoria(data))

    def test_toffoli_banderitas(self) -> None:
        r = theorem_k3_hash(b"abcdefgh" * 4, block_bits=32, banderitas=[0, 8])
        self.assertEqual(len(r.puntos_control), 2)
        self.assertLessEqual(r.final_hash, 0xFFFFFFFF)


if __name__ == "__main__":
    unittest.main()