#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Resolution.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_RESOLUTION_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_RESOLUTION_H_H

namespace antipc {

struct Resolution
{
    bool found = false;

    ResolutionSource source = ResolutionSource::None;

    ReferenceID reference = 0;

    float confidence = 0.0f;

    uint32_t estimatedCost = 0;
};

} // namespace antipc

#endif
