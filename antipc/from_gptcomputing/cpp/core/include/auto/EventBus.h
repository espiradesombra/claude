#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/EventBus.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_EVENTBUS_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_EVENTBUS_H_H

namespace antipc {

class EventBus
{
public:

    void push(Event e);

    bool empty() const;

    Event pop();

private:

    std::queue<Event> queue;
};

} // namespace antipc

#endif
