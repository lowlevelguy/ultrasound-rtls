#pragma once
#include "FS.h"
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <stdint.h>

#define LOG_FILE_PATH "/log.csv"
#define QUEUE_SIZE 10
#define START_BYTE 0xAA

extern SemaphoreHandle_t buffer_mutex;
extern TaskHandle_t logger_handle;

typedef struct logQueueEntry {
    float x;
    float y;
    unsigned long timestamp;
} log_queue_entry_t;

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

int initialize_log_file(fs::FS *fs);
void uart_recieve(void* param);
void logger_task(void* file_s);
