#pragma once
#include "FS.h"

#define LOG_FILE_PATH "/log.csv"
#define QUEUE_SIZE 10

typedef struct logQueue {
    float x;
    float y;
    unsigned long timestamp;
} log_queue_t;

int initialize_log_file(fs::FS &fs);
int save_to_queue(const float* pos, const unsigned long timestamp);
int save_to_file(fs::FS &fs);