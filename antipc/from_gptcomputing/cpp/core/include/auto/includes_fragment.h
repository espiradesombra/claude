// Fragmento includes lineas 1911-1918
#include "PayloadStore.h"

namespace antipc
{

PayloadID PayloadStore::store(std::vector<uint8_t>&& data)
{
    PayloadID id = m_nextId++;
