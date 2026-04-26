#pragma once

#include <stdint.h>
#include <Arduino.h>
#define MAX_WAIT_TIME 100

HardwareSerial *SENSORS[4];

int define_uart_port(int index);
int init_sensor(int index, int baud_rate);
uint16_t read_sensor(int sensor_index);