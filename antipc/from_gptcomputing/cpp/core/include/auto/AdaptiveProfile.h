#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/AdaptiveProfile.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_ADAPTIVEPROFILE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_ADAPTIVEPROFILE_H_H

namespace antipc {

struct AdaptiveProfile
{
    double averageExecutionTime;

    double averageReuse;

    double averageTrust;

    double averageMemory;

    double averageNetwork;

    double failureRate;
};

} // namespace antipc

#endif
