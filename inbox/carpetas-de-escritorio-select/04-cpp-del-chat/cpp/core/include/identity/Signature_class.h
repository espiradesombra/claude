#pragma once

/*
    AntiPC - Adaptive Network Through Parallel Computing
    ----------------------------------------------------

    Signature

    Representa la identidad criptográfica o lógica de un
    KnowledgeObject.

    El Runtime nunca asume cómo se genera.

    Puede representar:

      - SHA256
      - BLAKE3
      - XXHash
      - CRC
      - Firma experimental
      - Firma distribuida

    El algoritmo pertenece al Plugin.

    El Core solamente compara.
*/

#include <array>
#include <cstdint>

namespace antipc
{

class Signature
{
public:

    static constexpr std::size_t SIZE = 32;

    using Storage = std::array<std::uint8_t, SIZE>;

    constexpr Signature() = default;

    constexpr explicit Signature(const Storage& value)
        : data_(value)
    {
    }

    [[nodiscard]]
    constexpr const Storage& bytes() const noexcept
    {
        return data_;
    }

    [[nodiscard]]
    constexpr bool empty() const noexcept
    {
        for (auto b : data_)
        {
            if (b != 0)
                return false;
        }

        return true;
    }

    constexpr bool operator==(const Signature&) const = default;

private:

    Storage data_{};
};
