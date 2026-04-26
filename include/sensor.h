#pragma once

#include <stdint.h>
#include <Arduino.h>
#define MAX_WAIT_TIME 100
#define SENSOR_TRIGGER_BYTE 0x01

extern HardwareSerial *SENSORS[4];

int init_sensor(int sensor_index, int uart_index, int baud_rate, int RX, int TX);
uint16_t read_sensor(int sensor_index);