#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/IPayloadBackend.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_IPAYLOADBACKEND_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_IPAYLOADBACKEND_H_H

namespace antipc {

class IPayloadBackend
{
public:

    virtual PayloadID store(...) = 0;

    virtual bool load(...) = 0;

    virtual bool erase(...) = 0;

    virtual bool contains(...) = 0;

    virtual ~IPayloadBackend() = default;
};

} // namespace antipc

#endif
