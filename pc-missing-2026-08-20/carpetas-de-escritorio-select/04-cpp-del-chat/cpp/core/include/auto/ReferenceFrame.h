#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/ReferenceFrame.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_REFERENCEFRAME_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_REFERENCEFRAME_H_H

namespace antipc {

struct ReferenceFrame
{

ReferenceID reference;

PayloadID payload;

Signature signature;

Flags flags;

Tick tick;

}

} // namespace antipc

#endif
