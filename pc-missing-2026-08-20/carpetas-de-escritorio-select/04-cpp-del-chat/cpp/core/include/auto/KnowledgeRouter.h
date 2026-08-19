#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/KnowledgeRouter.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEROUTER_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_KNOWLEDGEROUTER_H_H

namespace antipc {

class KnowledgeRouter
{
public:

    Resolution resolve(const Signature& signature);

private:

    Resolution searchRam(const Signature&);

    Resolution searchKnowledge(const Signature&);

    Resolution searchPayload(const Signature&);

    Resolution searchShared(const Signature&);

    Resolution searchNetwork(const Signature&);

};

} // namespace antipc

#endif
