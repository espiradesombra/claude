#pragma once
// Extraido de gptcomputing.txt -> k3/include/KOPHeader.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_K3_INCLUDE_KOPHEADER_H_H
#define ANTIPC_K3_INCLUDE_KOPHEADER_H_H

namespace antipc {

Cada KOP tiene la misma estructura
struct KOPHeader
{
    uint16_t id;          // Tipo de operación
    uint16_t version;     // Versión del KOP
    uint32_t size;        // Tamaño
    uint64_t flags;       // Capacidades
};

} // namespace antipc

#endif

