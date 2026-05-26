#pragma once
#include <stdint.h>

typedef struct {
    uint8_t start;
    uint32_t timestamp;
    float pos[2];
    uint8_t checksum;
} __attribute__((packed)) log_packet_t;

typedef struct {
    uint8_t start;
    uint32_t timestamp;
} __attribute__((packed)) ack_packet_t;

int log_reading(uint32_t timestamp, float pos[2]);