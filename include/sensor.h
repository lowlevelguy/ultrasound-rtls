#pragma once

#include <stdint.h>
#include <Arduino.h>
#define MAX_WAIT_TIME_MS 100
#define SENSOR_TRIGGER_BYTE 0x01
#define MAX_SENSORS 4

extern HardwareSerial *SENSORS[MAX_SENSORS];

int init_sensor(int sensor_index, int uart_index, int baud_rate, int RX, int TX);
int32_t read_sensor(int sensor_index);