#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/KOPHeader.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_KOPHEADER_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_KOPHEADER_H_H

namespace antipc {

struct KOPHeader
{
    uint16_t id;          // Tipo de operación
    uint16_t version;     // Versión del KOP
    uint32_t size;        // Tamaño
    uint64_t flags;       // Capacidades
};

} // namespace antipc

#endif
