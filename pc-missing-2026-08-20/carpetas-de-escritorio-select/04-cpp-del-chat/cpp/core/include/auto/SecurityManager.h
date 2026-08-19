#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/SecurityManager.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_SECURITYMANAGER_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_SECURITYMANAGER_H_H

namespace antipc {

class SecurityManager
{
public:

    bool canRead(NodeID, ReferenceID);

    bool canWrite(NodeID, ReferenceID);

    bool canExecute(NodeID, PluginID);

    bool canShare(NodeID, ReferenceID);

};

} // namespace antipc

#endif
