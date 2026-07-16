#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Runtime.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_RUNTIME_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_RUNTIME_H_H

namespace antipc {

class Runtime
{
public:

    void start();

    void stop();

    void processNextEvent();

private:

    EventBus eventBus;

    Scheduler scheduler;

    PluginManager plugins;

    KnowledgeBuffer knowledge;

    DecisionEngine decision;
};

} // namespace antipc

#endif
