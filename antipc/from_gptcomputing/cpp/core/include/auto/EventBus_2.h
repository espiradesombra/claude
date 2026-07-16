#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/EventBus_2.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_EVENTBUS_2_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_EVENTBUS_2_H_H

namespace antipc {

class EventBus
{
public:

    bool push(const ReferenceFrame&);

    bool pop(ReferenceFrame&);

    std::size_t size() const;

private:

    RingBuffer<ReferenceFrame> queue;
};

} // namespace antipc

#endif
