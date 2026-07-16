#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/KnowledgeObject.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEOBJECT_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEOBJECT_H_H

namespace antipc {

class KnowledgeObject
{
Reference reference;

Payload payload;

Metadata metadata;

Statistics statistics;

History history;
};

} // namespace antipc

#endif
