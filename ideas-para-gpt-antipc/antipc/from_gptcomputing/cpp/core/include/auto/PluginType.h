#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/PluginType.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_PLUGINTYPE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_PLUGINTYPE_H_H

namespace antipc {

enum class PluginType
{
    Acquire,

    Transform,

    Infer,

    Verify,

    Publish,

    Effect
};

} // namespace antipc

#endif
