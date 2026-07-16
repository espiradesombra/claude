#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/PluginContract.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_PLUGINCONTRACT_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_PLUGINCONTRACT_H_H

namespace antipc {

struct PluginContract
{
    deterministic;

    idempotent;

    cacheable;

    distributable;

    sideEffects;

    estimatedCost;
};

} // namespace antipc

#endif
