#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/KOP.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_KOP_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_KOP_H_H

namespace antipc {

class KOP
{
public:
    virtual bool execute(Context&) = 0;
};

} // namespace antipc

#endif
