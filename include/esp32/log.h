#pragma once
#include "FS.h"
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <freertos/queue.h>
#include <freertos/semphr.h>
#include <LittleFS.h>
#include <stdint.h>


#define LOG_FILE_PATH "/history.log"
#define QUEUE_SIZE 10
#define START_BYTE 0xAA

extern TaskHandle_t logger_handle;
extern TaskHandle_t uart_handle;
extern uint32_t last_db_entry;
extern TaskHandle_t flush_to_db_handle;
extern SemaphoreHandle_t file_mutex;

typedef struct {
    uint32_t timestamp;
    float pos[2];
}timepos_t;

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


void uart_receive(void* param);
void logger_task(void* file_s);
int  cloud_begin();
int  cloud_publish(timepos_t* tp, size_t tp_size);
void cloud_loop();
void get_db_last_entry(void* param);
void flush_to_db(void* param);
void callback(char* topic, byte* message, unsigned int length);