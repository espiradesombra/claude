#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/PendingOperation.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_PENDINGOPERATION_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_PENDINGOPERATION_H_H

namespace antipc {

struct PendingOperation
{
    uint32_t required;

    uint32_t available;

    PluginID plugin;

    std::vector<ReferenceID> inputs;
};

} // namespace antipc

#endif
