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
inline float s1_correct(float meas) {
    return (meas + 21.5f) / 0.9896f;
}
inline float s2_correct(float meas) {
    return (meas + 12.6f) / 0.9876f;
}
inline float s3_correct(float meas) {
    return (meas + 19.9f) / 0.9713f;
}
inline float s4_correct(float meas) {
    return meas;
}
inline sensor_t sensors[4] = {
    {{.type = SerialHandle_t::HW, .hw = &Serial1}, 0x01, &s1_correct},
    {{.type = SerialHandle_t::HW, .hw = &Serial2}, 0x01, &s2_correct},
    {{.type = SerialHandle_t::HW, .hw = &Serial3}, 0x01, &s3_correct},
    {{.type = SerialHandle_t::HW, .hw = &Serial4}, 0x55, &s4_correct}
};

#define CONFIG_TLOG_SERIAL
inline HardwareSerial *tlog_serial = &Serial5;