#pragma once
#include "FS.h"
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>

#define LOG_FILE_PATH "/log.csv"
#define QUEUE_SIZE 10
TaskHandle_t logger_handle = NULL;

typedef struct logQueueEntry {
    float x;
    float y;
    unsigned long timestamp;
} log_queue_entry_t;

int initialize_log_file(fs::FS *fs);
int on_read(const float* pos, const unsigned long timestamp);
int save_to_file(fs::FS *fs, uint32_t buffer_index);