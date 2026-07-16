#pragma once
// Extraido de gptcomputing.txt -> core/src/knowledge/KnowledgeRouter_resolve.cpp
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_SRC_KNOWLEDGE_KNOWLEDGEROUTER_RESOLVE_CPP_H
#define ANTIPC_CORE_SRC_KNOWLEDGE_KNOWLEDGEROUTER_RESOLVE_CPP_H

namespace antipc {

Resolution searchPayload(const Signature&);

    Resolution searchShared(const Signature&);

    Resolution searchNetwork(const Signature&);

};

}
KnowledgeRouter.cpp

Lo primero será muy sencillo.

Resolution KnowledgeRouter::resolve(const Signature& s)
{
    Resolution r;

    r = searchRam(s);

    if(r.found)
        return r;

    r = searchKnowledge(s);

    if(r.found)
        return r;

    r = searchPayload(s);

    if(r.found)
        return r;

    r = searchShared(s);

    if(r.found)
        return r;

    r = searchNetwork(s);

    if(r.found)
        return r;

    r.source = ResolutionSource::ExecutePlugin;

    return r;

} // namespace antipc

#endif

