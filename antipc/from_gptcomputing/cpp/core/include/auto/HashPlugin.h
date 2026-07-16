#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/HashPlugin.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_HASHPLUGIN_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_HASHPLUGIN_H_H

namespace antipc {

class HashPlugin : public IPlugin
{
public:

    PluginInfo info() override;

    Signature signature() override;

    Cost estimate() override;

    ReferenceFrame execute(
        const ReferenceFrame& input
    ) override;
};

} // namespace antipc

#endif
