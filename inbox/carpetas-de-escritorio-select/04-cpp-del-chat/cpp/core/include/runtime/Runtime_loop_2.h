#pragma once
// Extraido de gptcomputing.txt -> core/include/runtime/Runtime_loop_2.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_RUNTIME_RUNTIME_LOOP_2_H_H
#define ANTIPC_CORE_INCLUDE_RUNTIME_RUNTIME_LOOP_2_H_H

namespace antipc {

while(runtime.running())
{
    runtime.processNextEvent();
}

} // namespace antipc

#endif
