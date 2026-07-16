#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/ResolutionSource.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_RESOLUTIONSOURCE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_RESOLUTIONSOURCE_H_H

namespace antipc {

enum class ResolutionSource
{
    None,

    RAM,

    KnowledgeBuffer,

    PayloadStore,

    SharedMemory,

    Network,

    ExecutePlugin
};

} // namespace antipc

#endif
