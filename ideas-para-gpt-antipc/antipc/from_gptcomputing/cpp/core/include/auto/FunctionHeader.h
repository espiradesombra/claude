#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/FunctionHeader.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_FUNCTIONHEADER_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_FUNCTIONHEADER_H_H

namespace antipc {

struct FunctionHeader
{
    uint16_t id;
    uint16_t version;
    uint32_t size;
};

} // namespace antipc

#endif
