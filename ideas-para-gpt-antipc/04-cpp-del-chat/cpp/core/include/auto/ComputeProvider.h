#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/ComputeProvider.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_COMPUTEPROVIDER_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_COMPUTEPROVIDER_H_H

namespace antipc {

class ComputeProvider
{
public:

    virtual Capability capability();

    virtual Cost estimate();

    virtual Result execute();
};

} // namespace antipc

#endif
