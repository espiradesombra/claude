#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Capability_2.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_CAPABILITY_2_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_CAPABILITY_2_H_H

namespace antipc {

struct Capability
{
    bool deterministic;

    bool parallel;

    bool distributable;

    bool cacheable;

    bool requires_network;

    bool requires_permission;
};

} // namespace antipc

#endif
