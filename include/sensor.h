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

extern SerialHandle_t SENSORS[4];

int init_sensor(int sensor_index, SerialHandle_t* handle, int baud_rate, int RX, int TX);
int32_t read_sensor(int sensor_index);