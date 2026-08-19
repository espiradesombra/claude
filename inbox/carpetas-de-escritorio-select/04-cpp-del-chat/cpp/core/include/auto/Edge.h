#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Edge.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_EDGE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_EDGE_H_H

namespace antipc {

struct Edge
{
    ReferenceID from;
    ReferenceID to;
    PluginID plugin;
};

} // namespace antipc

#endif
