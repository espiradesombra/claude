#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/KnowledgeDNA_2.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEDNA_2_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEDNA_2_H_H

namespace antipc {

struct KnowledgeDNA
{
    // Identidad lógica del contenido.
    Signature signature;

    // Plugin productor.
    PluginID producer = 0;

    // Versión del plugin.
    Version producerVersion = 0;

    // Política utilizada.
    PolicyID policy = 0;

    // Conocimiento origen.
    KnowledgeID parent = INVALID_KNOWLEDGE_ID;

    // Época lógica.
    std::uint64_t epoch = 0;

    // Tick interno.
    std::uint64_t tick = 0;

    [[nodiscard]]
    constexpr bool hasParent() const noexcept
    {
        return isValid(parent);
    }
};

} // namespace antipc

#endif
