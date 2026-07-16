#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/ReferenceState.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_REFERENCESTATE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_REFERENCESTATE_H_H

namespace antipc {

enum class ReferenceState : uint8_t
{
    Created,
    Unverified,
    Verified,
    Published,
    Active,
    Archived,
    Invalid
};

} // namespace antipc

#endif
