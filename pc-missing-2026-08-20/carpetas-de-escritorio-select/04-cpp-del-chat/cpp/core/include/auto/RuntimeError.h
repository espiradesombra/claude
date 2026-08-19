#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/RuntimeError.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_RUNTIMEERROR_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_RUNTIMEERROR_H_H

namespace antipc {

enum class RuntimeError
{
    None,

    PluginNotFound,

    InvalidReference,

    InvalidPayload,

    PermissionDenied,

    NetworkUnavailable,

    VerificationFailed
};

} // namespace antipc

#endif
