#pragma once
// Extraido de gptcomputing.txt -> core/src/storage/PayloadStore.cpp
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_SRC_STORAGE_PAYLOADSTORE_CPP_H
#define ANTIPC_CORE_SRC_STORAGE_PAYLOADSTORE_CPP_H

namespace antipc {

#include "PayloadStore.h"

namespace antipc
{

PayloadID PayloadStore::store(std::vector<uint8_t>&& data)
{
    PayloadID id = m_nextId++;

    m_storage.emplace(id,std::move(data));

    return id;
}

const std::vector<uint8_t>* PayloadStore::get(PayloadID id) const
{
    auto it = m_storage.find(id);

    if(it==m_storage.end())
        return nullptr;

    return &it->second;
}

bool PayloadStore::exists(PayloadID id) const
{
    return m_storage.contains(id);
}

void PayloadStore::remove(PayloadID id)
{
    m_storage.erase(id);
}

std::size_t PayloadStore::size() const
{
    return m_storage.size();
}

}

} // namespace antipc

#endif

