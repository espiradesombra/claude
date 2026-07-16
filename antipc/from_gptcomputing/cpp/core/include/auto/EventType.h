#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/EventType.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_EVENTTYPE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_EVENTTYPE_H_H

namespace antipc {

enum class EventType
{
    CreateReference,

    ExecuteOperation,

    Publish,

    Verify,

    Shutdown
};

} // namespace antipc

#endif
