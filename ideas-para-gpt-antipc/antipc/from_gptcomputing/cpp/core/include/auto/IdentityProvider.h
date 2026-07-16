#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/IdentityProvider.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_IDENTITYPROVIDER_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_IDENTITYPROVIDER_H_H

namespace antipc {

class IdentityProvider
{
public:

    virtual K3Identity create(...) = 0;

    virtual bool verify(...) = 0;
};

} // namespace antipc

#endif
