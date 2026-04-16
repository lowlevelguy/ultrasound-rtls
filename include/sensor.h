#pragma once

#include <stdint.h>

extern const int SENSORS[4];

uint64_t read_sensor(int sensor_index);