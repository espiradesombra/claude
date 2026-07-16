#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/PluginInfo.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_PLUGININFO_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_PLUGININFO_H_H

namespace antipc {

struct PluginInfo
{
    std::string name;

    std::string author;

    uint32_t version;

    uint32_t api_version;

    PluginType type;
};

} // namespace antipc

#endif
