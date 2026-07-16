#pragma once

#include "Signature.h"
#include "Reference.h"

namespace antipc
{

enum class ResolutionSource
{
    None,

    RAM,

    KnowledgeBuffer,

    PayloadStore,

    SharedMemory,

    Network,

    ExecutePlugin
};

struct Resolution
{
    bool found = false;

    ResolutionSource source = ResolutionSource::None;

    ReferenceID reference = 0;

    float confidence = 0.0f;
