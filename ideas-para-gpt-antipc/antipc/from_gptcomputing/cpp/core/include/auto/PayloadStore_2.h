#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/PayloadStore_2.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_PAYLOADSTORE_2_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_PAYLOADSTORE_2_H_H

namespace antipc {

class PayloadStore
{
public:

    PayloadID store(std::vector<uint8_t>&& data);

    const std::vector<uint8_t>* get(PayloadID id) const;

    bool exists(PayloadID id) const;

    void remove(PayloadID id);

    std::size_t size() const;

private:

    PayloadID m_nextId = 1;

    std::unordered_map<PayloadID,
                       std::vector<uint8_t>> m_storage;
};

} // namespace antipc

#endif
