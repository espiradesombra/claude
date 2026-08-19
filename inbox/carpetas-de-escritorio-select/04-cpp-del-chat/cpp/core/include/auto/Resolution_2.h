#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Resolution_2.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_RESOLUTION_2_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_RESOLUTION_2_H_H

namespace antipc {

enum class Resolution
{
    AcceptHighestTrust,

    RequestVerification,

    ExecuteAgain,

    Wait,

    Reject,

    ManualDecision
};

} // namespace antipc

#endif
