#pragma once

#include <stdint.h>

int  cloud_begin();
int  cloud_publish(uint32_t timestamp, const float pos[2]);
void cloud_loop();