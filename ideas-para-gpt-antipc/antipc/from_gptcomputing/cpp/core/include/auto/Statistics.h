#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Statistics.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_STATISTICS_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_STATISTICS_H_H

namespace antipc {

struct Statistics
{
    uint64_t hits = 0;

    Tick created = 0;

    Tick verified = 0;

    Tick last_access = 0;
};

} // namespace antipc

#endif
