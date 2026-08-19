#pragma once

/*
    AntiPC - Adaptive Network Through Parallel Computing
    ----------------------------------------------------

    KnowledgeDNA

    Describe el origen de un KnowledgeObject.

    El contenido puede ser idéntico entre dos objetos,
    pero su ADN puede ser distinto si fueron obtenidos
    mediante procesos diferentes.

    El DNA nunca cambia.
*/

#include "KnowledgeID.h"
#include "Signature.h"

#include <cstdint>

namespace antipc
{

using PluginID = std::uint32_t;
using PolicyID = std::uint32_t;
using Version  = std::uint32_t;

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

}
