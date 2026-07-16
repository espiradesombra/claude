#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/KnowledgeValue.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEVALUE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEVALUE_H_H

namespace antipc {

struct KnowledgeValue
{
    double creationCost;
    double storageCost;
    double transferCost;

    uint64_t reuseCount;

    uint64_t lastUse;

    double trust;

    double probabilityReuse;

    double size;
};

} // namespace antipc

#endif
