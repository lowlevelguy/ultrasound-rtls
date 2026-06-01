#pragma once
#include <stdint.h>

#define START_BYTE 0xAA

#define ERROR_LOG_WRITE_TIMEOUT     -1
#define ERROR_LOG_WRITE_FAILURE     -2
#define ERROR_LOG_ACK_TIMEOUT       -3
#define ERROR_LOG_READ_FAILURE      -4
#define ERROR_LOG_ACK_INVALID       -5

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

int log_timepos(uint32_t timestamp, float pos[2]);