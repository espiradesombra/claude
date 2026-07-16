#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Metadata.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_METADATA_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_METADATA_H_H

namespace antipc {

struct Metadata
{
    float confidence = 0.0f;

    PluginID creator = 0;

    uint32_t version = 1;

    uint32_t flags = 0;
};

} // namespace antipc

#endif
