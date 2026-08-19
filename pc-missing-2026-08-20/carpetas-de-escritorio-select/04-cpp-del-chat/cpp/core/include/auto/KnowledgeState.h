#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/KnowledgeState.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGESTATE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGESTATE_H_H

namespace antipc {

enum class KnowledgeState
{
    Created,

    Verified,

    Cached,

    Shared,

    Archived,

    Invalid
};

} // namespace antipc

#endif
