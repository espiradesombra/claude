/**
 * ZypyZape / AntiPC — Maestro UDP Windows (viabilidad)
 * Víctor Manzanares Alberola · VMA 2026
 *
 * HONESTO: recvfrom(buf) NO es DMA directo NIC→RAM arbitraria.
 * El kernel copia driver→espacio usuario UNA vez. Lo que evitamos es la
 * SEGUNDA copia (temporal→cola) usando slots fijos del ring.
 */
#define WIN32_LEAN_AND_MEAN
#include <winsock2.h>
#include <ws2tcpip.h>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <thread>
#include <vector>

#pragma comment(lib, "ws2_32.lib")

static constexpr size_t BUFFER_SIZE = 1024;
static constexpr size_t DATA_SIZE   = 64;
static constexpr int    PORT        = 3333;

struct alignas(64) Packet {
    uint8_t payload[DATA_SIZE];
};

static std::vector<Packet> g_ring(BUFFER_SIZE);
static std::atomic<size_t> g_head{0};
static std::atomic<size_t> g_tail{0};
static std::atomic<bool>   g_running{true};
static std::atomic<uint64_t> g_received{0};
static std::atomic<uint64_t> g_processed{0};
static std::atomic<uint64_t> g_drops{0};

static void network_ingest() {
    WSADATA wsa{};
    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0) return;

    SOCKET sock = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    u_long nb = 1;
    ioctlsocket(sock, FIONBIO, &nb);

    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_port = htons(PORT);
    addr.sin_addr.s_addr = INADDR_ANY;
    bind(sock, reinterpret_cast<sockaddr*>(&addr), sizeof(addr));

    std::cout << "[RED] UDP :" << PORT << " non-blocking, recvinto ring slots\n";

    char discard[2048];
    while (g_running.load(std::memory_order_relaxed)) {
        size_t h = g_head.load(std::memory_order_relaxed);
        size_t next = (h + 1) % BUFFER_SIZE;
        if (next == g_tail.load(std::memory_order_acquire)) {
            sockaddr_in from{};
            int flen = sizeof(from);
            recvfrom(sock, discard, sizeof(discard), 0,
                     reinterpret_cast<sockaddr*>(&from), &flen);
            g_drops.fetch_add(1, std::memory_order_relaxed);
            std::this_thread::sleep_for(std::chrono::microseconds(50));
            continue;
        }

        sockaddr_in from{};
        int flen = sizeof(from);
        int n = recvfrom(sock,
                         reinterpret_cast<char*>(g_ring[h].payload),
                         static_cast<int>(DATA_SIZE),
                         0,
                         reinterpret_cast<sockaddr*>(&from),
                         &flen);
        if (n > 0) {
            g_head.store(next, std::memory_order_release);
            g_received.fetch_add(1, std::memory_order_relaxed);
        } else {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }
    closesocket(sock);
    WSACleanup();
}

static void processing_channel() {
    while (g_running.load(std::memory_order_relaxed)) {
        size_t t = g_tail.load(std::memory_order_relaxed);
        if (t == g_head.load(std::memory_order_acquire)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            continue;
        }
        const Packet& pkt = g_ring[t];
        uint32_t acc = 0;
        for (size_t i = 0; i < DATA_SIZE; ++i)
            acc = (acc * 131u) ^ pkt.payload[i];
        (void)acc;
        g_tail.store((t + 1) % BUFFER_SIZE, std::memory_order_release);
        g_processed.fetch_add(1, std::memory_order_relaxed);
    }
}

int main() {
    std::cout << "ZypyZape master — slot ring UDP\n";
    std::thread rx(network_ingest);
    std::thread logic(processing_channel);

    std::this_thread::sleep_for(std::chrono::seconds(10));
    g_running = false;
    rx.join();
    logic.join();

    std::cout << "recv=" << g_received.load()
              << " proc=" << g_processed.load()
              << " drops=" << g_drops.load() << "\n";
    return 0;
}