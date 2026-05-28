#pragma once
#include "teensy/sensor.h"

#define CONFIG_ANCHOR_POS
inline constexpr float anchor_pos[4][2] = {
    {0, 0},
    {4000, 0},
    {0, 4000},
    {4000, 4000}
};

#define CONFIG_SENSORS
inline SerialHandle_t SENSORS[4] = {
    {.type = SerialHandle_t::HW, .hw = &Serial1},
    {.type = SerialHandle_t::HW, .hw = &Serial2},
    {.type = SerialHandle_t::HW, .hw = &Serial3},
    {.type = SerialHandle_t::HW, .hw = &Serial4}
};

#define CONFIG_TLOG_SERIAL
inline HardwareSerial *tlog_serial = &Serial5;