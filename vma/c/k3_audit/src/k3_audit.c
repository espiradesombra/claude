/*
 * k3_audit.c — Motor K3 industrial con firma activa, stride y marcas finales.
 */

#include "k3_audit.h"

#include <string.h>

const char K3_FIRMA_LOGICA[] = "x^=(B*0x9E3779B1);x+=rotar(B,6);estado[0]=x;";

static const int OFFSETS[8] = {5, 7, 9, 11, 13, 17, 19, 23};

static inline uint32_t rotar_izquierda_32bits(uint32_t valor, int posiciones) {
    posiciones &= 31;
    if (posiciones == 0) {
        return valor;
    }
    return (valor << posiciones) | (valor >> (32 - posiciones));
}

static void ejecutar_compresion_k3(uint32_t* estado, uint32_t valor_bloque,
                                   const K3AuditConfig* config) {
    int N = config->num_registros;
    uint32_t magico = config->constante_disp;
    uint32_t B = valor_bloque;

    uint32_t x = estado[0] ^ (B & estado[1 % N]);
    int j;
    for (j = 2; j < N; j++) {
        x ^= rotar_izquierda_32bits(estado[j], OFFSETS[j % 8]);
    }
    x ^= rotar_izquierda_32bits(estado[1 % N], OFFSETS[4]);
    x += rotar_izquierda_32bits(estado[0], OFFSETS[5]);
    x ^= (B * magico);

    x ^= B;
    x += rotar_izquierda_32bits(B, OFFSETS[6]);
    x ^= rotar_izquierda_32bits(B, OFFSETS[7]);
    x ^= rotar_izquierda_32bits(x, OFFSETS[2]);
    x *= magico;
    x ^= (x >> 15);

    uint32_t estado_anterior_0 = estado[0];
    estado[0] = x;

    uint32_t temp_prev = estado_anterior_0;
    int i;
    for (i = 1; i < N; i++) {
        uint32_t temp_actual = estado[i];
        estado[i] = rotar_izquierda_32bits(estado[i], OFFSETS[0] + i) ^ temp_prev;
        temp_prev = temp_actual;
        if (i == N - 1) {
            estado[i] ^= B;
        }
    }
}

static uint32_t extraer_hash_final(const uint32_t* estado_final, const K3AuditConfig* config) {
    uint32_t acumulador = estado_final[0];
    int i;
    for (i = 1; i < config->num_registros; i++) {
        acumulador ^= rotar_izquierda_32bits(estado_final[i], 5 + i);
    }
    acumulador ^= acumulador >> 16;
    acumulador *= 0x85EBCA6Bu;
    acumulador ^= acumulador >> 13;
    acumulador *= 0xC2B2AE35u;
    acumulador ^= acumulador >> 16;
    return acumulador;
}

K3AuditConfig k3_audit_config_default(void) {
    return k3_audit_crear_config(32, 0, 3, 0x1F2E3D4Cu, 0x11223344u);
}

K3AuditConfig k3_audit_crear_config(int bits_lectura, int desfase_trozos,
                                    int n_registros, uint32_t semilla,
                                    uint32_t marca_inicial) {
    K3AuditConfig config;
    int i;

    if (bits_lectura != 8 && bits_lectura != 16 && bits_lectura != 32) {
        config.bits_ancho = 32;
    } else {
        config.bits_ancho = bits_lectura;
    }

    config.bytes_ancho = config.bits_ancho / 8;
    config.desfase_stride = (desfase_trozos < 0) ? 0 : desfase_trozos;
    config.num_registros = (n_registros < 2) ? 2 : n_registros;
    if (config.num_registros > K3_AUDIT_MAX_REGISTROS) {
        config.num_registros = K3_AUDIT_MAX_REGISTROS;
    }
    config.semilla_inicial = semilla;
    config.marca_inicial = marca_inicial;
    config.constante_disp = 0x9E3779B1u;
    config.total_marcas_encoladas = 0;

    for (i = 0; i < K3_AUDIT_MAX_MARCAS; i++) {
        config.cola_marcas[i].activa = 0;
    }

    return config;
}

void k3_audit_encolar_marca(K3AuditConfig* config, uint32_t valor_bits, uint8_t banderita) {
    int idx;
    if (!config) {
        return;
    }
    idx = config->total_marcas_encoladas;
    if (idx >= K3_AUDIT_MAX_MARCAS) {
        return;
    }
    config->cola_marcas[idx].valor_bits = valor_bits;
    config->cola_marcas[idx].banderita = banderita;
    config->cola_marcas[idx].activa = 1;
    config->total_marcas_encoladas++;
}

void k3_audit_encolar_marcas_default(K3AuditConfig* config) {
    k3_audit_encolar_marca(config, 0x00010203u, 0x10);
    k3_audit_encolar_marca(config, 0x55555555u, 0xEE);
}

uint32_t k3_audit_hash_buffer(const uint8_t* datos, size_t longitud,
                              const K3AuditConfig* config,
                              uint8_t* telemetria_out, size_t telemetria_cap,
                              size_t* puntos_telemetria,
                              uint32_t* estado_out) {
    int N;
    size_t contador_telemetria = 0;
    int ancho_bytes;
    int desfase;
    size_t len_firma;
    size_t idx_firma;
    size_t indice;
    int i;
    uint32_t estado_local[K3_AUDIT_MAX_REGISTROS];

    if (!config || !datos) {
        return 0;
    }

    N = config->num_registros;
    ancho_bytes = config->bytes_ancho;
    desfase = config->desfase_stride;

    for (i = 0; i < N; i++) {
        estado_local[i] = config->semilla_inicial ^ (i * config->constante_disp);
    }

    len_firma = strlen(K3_FIRMA_LOGICA);
    idx_firma = 0;
    while (idx_firma < len_firma) {
        uint32_t bloque_firma = 0;
        int bytes_leidos = 0;
        int b;
        for (b = 0; b < ancho_bytes; b++) {
            if (idx_firma + (size_t)b < len_firma) {
                bloque_firma = (bloque_firma << 8) | (uint8_t)K3_FIRMA_LOGICA[idx_firma + b];
                bytes_leidos++;
            }
        }
        if (bytes_leidos == 0) {
            break;
        }
        if (bytes_leidos < ancho_bytes) {
            bloque_firma <<= (8 * (ancho_bytes - bytes_leidos));
        }
        ejecutar_compresion_k3(estado_local, bloque_firma, config);
        if (telemetria_out && contador_telemetria < telemetria_cap) {
            uint32_t paridad = 0;
            for (i = 0; i < N; i++) {
                paridad ^= estado_local[i];
            }
            telemetria_out[contador_telemetria++] = (uint8_t)(paridad & 0xFF);
        }
        idx_firma += (size_t)bytes_leidos;
    }

    indice = 0;
    while (indice < longitud) {
        uint32_t valor_bloque = 0;
        int bytes_leidos = 0;
        int b;
        for (b = 0; b < ancho_bytes; b++) {
            if (indice + (size_t)b < longitud) {
                valor_bloque = (valor_bloque << 8) | datos[indice + b];
                bytes_leidos++;
            }
        }
        if (bytes_leidos == 0) {
            break;
        }
        if (bytes_leidos < ancho_bytes) {
            valor_bloque <<= (8 * (ancho_bytes - bytes_leidos));
        }
        ejecutar_compresion_k3(estado_local, valor_bloque, config);
        if (telemetria_out && contador_telemetria < telemetria_cap) {
            uint32_t paridad = 0;
            for (i = 0; i < N; i++) {
                paridad ^= estado_local[i];
            }
            telemetria_out[contador_telemetria++] =
                (uint8_t)((paridad ^ valor_bloque) & 0xFF);
        }
        indice += (size_t)ancho_bytes + (size_t)desfase;
    }

    for (i = 0; i < K3_AUDIT_MAX_MARCAS; i++) {
        if (config->cola_marcas[i].activa) {
            uint32_t bloque_compuesto = config->cola_marcas[i].valor_bits ^
                ((uint32_t)config->cola_marcas[i].banderita << 24);
            ejecutar_compresion_k3(estado_local, bloque_compuesto, config);
            if (telemetria_out && contador_telemetria < telemetria_cap) {
                uint32_t paridad = 0;
                int r;
                for (r = 0; r < N; r++) {
                    paridad ^= estado_local[r];
                }
                telemetria_out[contador_telemetria++] =
                    (uint8_t)((paridad ^ bloque_compuesto) & 0xFF);
            }
        }
    }

    if (puntos_telemetria) {
        *puntos_telemetria = contador_telemetria;
    }
    if (estado_out) {
        for (i = 0; i < N; i++) {
            estado_out[i] = estado_local[i];
        }
    }

    return extraer_hash_final(estado_local, config);
}