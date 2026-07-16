#pragma once
// Extraido de gptcomputing.txt -> k3/include/KOP.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_K3_INCLUDE_KOP_H_H
#define ANTIPC_K3_INCLUDE_KOP_H_H

namespace antipc {

class KOP
{
public:
    virtual bool execute(Context&) = 0;
};

} // namespace antipc

#endif

