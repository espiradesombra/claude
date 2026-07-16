#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/SignatureRecord.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_SIGNATURERECORD_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_SIGNATURERECORD_H_H

namespace antipc {

struct SignatureRecord
{
    Signature digest;

    AlgorithmID algorithm;

    uint64_t timestamp;

    uint32_t revision;

    uint32_t signatureCount;
};

} // namespace antipc

#endif
