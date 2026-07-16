#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Worker.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_WORKER_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_WORKER_H_H

namespace antipc {

class Worker
{
public:

    void execute(const ReferenceFrame&);
};

} // namespace antipc

#endif
