#pragma once
// Extraido de gptcomputing.txt -> main.cpp
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_MAIN_CPP_H
#define ANTIPC_MAIN_CPP_H

namespace antipc {

int main()
{
    Runtime runtime;
    runtime.loadPlugin("Demo");
    runtime.start();
    runtime.publish(5);
    runtime.run();
    return 0;
}

} // namespace antipc

#endif
