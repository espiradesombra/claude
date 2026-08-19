#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/K3Identity.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_K3IDENTITY_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_K3IDENTITY_H_H

namespace antipc {

struct K3Identity
{
    Signature k3Signature;

    KnowledgeDNA dna;

    uint64_t logicalTime;

    uint32_t revision;

    uint32_t trust;

    uint32_t flags;

    std::vector<Signature> externalSignatures;
};

} // namespace antipc

#endif
