#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Cost.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_COST_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_COST_H_H

namespace antipc {

struct Cost
{
    double cpu;

    double latency;

    double bandwidth;

    double energy;

    double confidence;
};

} // namespace antipc

#endif
