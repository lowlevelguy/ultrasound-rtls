#pragma once

#include <stdint.h>
#include <Arduino.h>
#include <SoftwareSerial.h>
#include <HardwareSerial.h>
#include "freertos/FreeRTOS.h"
#include "freertos/task.h"

#define MAX_WAIT_TIME_MS 100
#define SENSOR_TRIGGER_BYTE 0x01
#define MAX_SENSORS 4

typedef struct SerialHandle {
    enum { HW, SW } type;
    union {
        HardwareSerial* hw;
        SoftwareSerial* sw;
    };
} SerialHandle_t;

extern SerialHandle_t *SENSORS[MAX_SENSORS];

int init_sensor(int sensor_index, SerialHandle_t* handle, int baud_rate, int RX, int TX);
int32_t read_sensor(int sensor_index);