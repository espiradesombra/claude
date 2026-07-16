#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Signature.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_SIGNATURE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_SIGNATURE_H_H

namespace antipc {

struct Signature
{
    std::array<uint8_t,32> value{};

    bool operator==(const Signature&) const = default;
};

} // namespace antipc

#endif
