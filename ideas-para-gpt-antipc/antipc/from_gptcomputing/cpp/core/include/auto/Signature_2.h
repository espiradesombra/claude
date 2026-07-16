#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Signature_2.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_SIGNATURE_2_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_SIGNATURE_2_H_H

namespace antipc {

struct Signature
{
    uint64_t high;
    uint64_t low;

    bool operator==(const Signature&) const = default;
};

} // namespace antipc

#endif
