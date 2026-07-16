#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/KnowledgeCost.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGECOST_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGECOST_H_H

namespace antipc {

struct KnowledgeCost
{
    double cpu;

    double memory;

    double storage;

    double bandwidth;

    double latency;

    double energy;
};

} // namespace antipc

#endif
