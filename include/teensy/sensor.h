#pragma once

#include <stdint.h>
#include <SoftwareSerial.h>
#include <HardwareSerial.h>

typedef struct SerialHandle {
    enum { HW, SW } type;
    union {
        HardwareSerial* hw;
        SoftwareSerial* sw;
    };
} SerialHandle_t;

int sensor_begin(int sensor_index, int baud_rate);
int32_t read_sensor(int sensor_index);