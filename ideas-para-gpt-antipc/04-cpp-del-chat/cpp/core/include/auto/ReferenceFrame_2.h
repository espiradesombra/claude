#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/ReferenceFrame_2.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_REFERENCEFRAME_2_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_REFERENCEFRAME_2_H_H

namespace antipc {

struct ReferenceFrame
{
    ReferenceID reference;

    KnowledgeID knowledge;

    PayloadID payload;

    Signature signature;

    PluginID producer;

    Tick tick;

    uint32_t flags;

    float confidence;
};

} // namespace antipc

#endif
