#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/KnowledgeDNA.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEDNA_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEDNA_H_H

namespace antipc {

struct KnowledgeDNA
{
    Signature signature;

    PluginID producer;

    Version pluginVersion;

    PolicyID policy;

    std::vector<ReferenceID> parents;
};

} // namespace antipc

#endif
