#pragma once

/*
    AntiPC - Adaptive Network Through Parallel Computing
    ----------------------------------------------------

    KnowledgeID

    Identificador único de un objeto de conocimiento.

    Reglas:
      - Nunca cambia una vez asignado.
      - Es válido en todo el Runtime.
      - No contiene significado semántico.
      - Puede mapearse posteriormente a UUID, hash u otro formato.

    Versión:
      v0.1.0-alpha
*/

#include <cstdint>
#include <limits>

namespace antipc
{

using KnowledgeID = std::uint64_t;

constexpr KnowledgeID INVALID_KNOWLEDGE_ID =
    std::numeric_limits<KnowledgeID>::max();

constexpr bool isValid(KnowledgeID id)
{
    return id != INVALID_KNOWLEDGE_ID;
}

} // namespace antipc
