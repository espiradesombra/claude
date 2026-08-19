#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/IPlugin.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_IPLUGIN_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_IPLUGIN_H_H

namespace antipc {

class IPlugin
{
public:

    virtual PluginInfo info() = 0;

    virtual Signature signature() = 0;

    virtual Cost estimate() = 0;

    virtual bool validate() = 0;

    virtual Reference execute() = 0;

    virtual ~IPlugin() = default;
};

} // namespace antipc

#endif
