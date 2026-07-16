#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/KnowledgeFrame.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEFRAME_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEFRAME_H_H

namespace antipc {

struct KnowledgeFrame
{
    Header;

    Payload;

    Evidence;

    Policy;

    Context;

    Metadata;
};

} // namespace antipc

#endif
