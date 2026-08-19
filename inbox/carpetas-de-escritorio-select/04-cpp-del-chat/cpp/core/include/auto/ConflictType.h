#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/ConflictType.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_CONFLICTTYPE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_CONFLICTTYPE_H_H

namespace antipc {

enum class ConflictType
{
    None,

    DifferentPayload,

    DifferentEvidence,

    DifferentTrust,

    DifferentPolicy,

    DifferentPlugin,

    NetworkDisagreement,

    Timeout,

    Unknown
};

} // namespace antipc

#endif
