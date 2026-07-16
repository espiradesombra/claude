#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Capability.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_CAPABILITY_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_CAPABILITY_H_H

namespace antipc {

struct Capability
{
    bool deterministic;

    bool distributable;

    bool cacheable;

    bool reversible;

    bool streamable;

    bool requires_network;

    bool requires_permission;
};

} // namespace antipc

#endif
