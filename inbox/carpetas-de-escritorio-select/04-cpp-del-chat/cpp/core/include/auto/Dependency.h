#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Dependency.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_DEPENDENCY_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_DEPENDENCY_H_H

namespace antipc {

struct Dependency
{
    PluginID id;

    uint32_t minimum_version;
};

} // namespace antipc

#endif
