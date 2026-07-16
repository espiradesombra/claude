#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/ValueEngine.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_VALUEENGINE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_VALUEENGINE_H_H

namespace antipc {

class ValueEngine
{
public:

    KnowledgeDecision evaluate(
        const KnowledgeObject&
    );
};

} // namespace antipc

#endif
