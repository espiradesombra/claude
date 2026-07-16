#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace antipc
{

using Tick = uint64_t;

using ReferenceID = uint64_t;
using KnowledgeID = uint64_t;
using PluginID = uint32_t;

struct Signature
{
    std::array<uint8_t,32> value{};

    bool operator==(const Signature&) const = default;
};

enum class ReferenceState : uint8_t
{
    Created,
    Unverified,
    Verified,
    Published,
    Active,
    Archived,
    Invalid
};

struct Statistics
{
    uint64_t hits = 0;

    Tick created = 0;

    Tick verified = 0;

    Tick last_access = 0;
};

struct Metadata
{
    float confidence = 0.0f;

    PluginID creator = 0;

    uint32_t version = 1;

    uint32_t flags = 0;
};

class Reference
{
public:

    Reference() = default;

    explicit Reference(ReferenceID id)
        :
        m_id(id)
    {}

    ReferenceID id() const
    {
        return m_id;
    }

    KnowledgeID knowledge() const
    {
        return m_knowledge;
    }

    const Signature& signature() const
    {
        return m_signature;
    }

    ReferenceState state() const
    {
        return m_state;
    }

    Metadata& metadata()
    {
        return m_metadata;
    }

    Statistics& statistics()
    {
        return m_stats;
    }

    const std::vector<ReferenceID>& parents() const
    {
        return m_parents;
    }

    const std::vector<ReferenceID>& children() const
    {
        return m_children;
    }

    void addParent(ReferenceID id)
    {
        m_parents.push_back(id);
    }

    void addChild(ReferenceID id)
    {
        m_children.push_back(id);
    }

    void setState(ReferenceState s)
    {
        m_state = s;
    }

    void touch(Tick tick)
    {
        ++m_stats.hits;
        m_stats.last_access = tick;
    }

private:

    ReferenceID m_id = 0;

    KnowledgeID m_knowledge = 0;

    Signature m_signature;

    ReferenceState m_state = ReferenceState::Created;

    Metadata m_metadata;

    Statistics m_stats;

    std::vector<ReferenceID> m_parents;

    std::vector<ReferenceID> m_children;
};

}
