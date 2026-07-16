#pragma once
// Extraido de gptcomputing.txt -> core/include/auto/Reference.h
// Implementacion Python activa: src/antipc/runtime/

#ifndef ANTIPC_CORE_INCLUDE_AUTO_REFERENCE_H_H
#define ANTIPC_CORE_INCLUDE_AUTO_REFERENCE_H_H

namespace antipc {

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

} // namespace antipc

#endif
