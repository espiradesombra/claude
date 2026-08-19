#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/RingBuffer.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_RINGBUFFER_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_RINGBUFFER_H_H

namespace antipc {

template<typename T,std::size_t N>

class RingBuffer
{
public:

    bool push(const T& value);

    bool pop(T& value);

private:

    std::array<T,N> buffer;

    std::atomic<uint32_t> head{0};

    std::atomic<uint32_t> tail{0};
};

} // namespace antipc

#endif
