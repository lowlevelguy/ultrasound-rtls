#include <math.h>
#include "teensy/sensor.h"
#include "teensy/config.h"

#ifndef CONFIG_SENSORS
#error "Please define log_serial in config.h"
#endif

#define MAX_WAIT_TIME_MS 100

static int sensor_available(SerialHandle_t s) {
    return (s.type == SerialHandle_t::HW)
        ? s.hw->available()
        : s.sw->available();
}

static size_t sensor_read(SerialHandle_t s) {
    return (s.type == SerialHandle_t::HW)
        ? s.hw->read()
        : s.sw->read();
}

static void sensor_write(SerialHandle_t s, uint8_t data) {
    if (s.type == SerialHandle_t::HW)
        s.hw->write(data);
    else
        s.sw->write(data);
}

static size_t sensor_read_bytes(SerialHandle_t s, uint8_t* buffer, size_t len) {
    if (s.type == SerialHandle_t::HW)
        return s.hw->readBytes(buffer, len);
    else
        return s.sw->readBytes(buffer, len);
}

int sensor_begin(int sensor_index, int baud_rate) {
    if (sensor_index < 0 || sensor_index >= 4)
        return -1;

    if (sensors[sensor_index].serial.type == SerialHandle_t::HW)
        sensors[sensor_index].serial.hw->begin(baud_rate);
    else
        sensors[sensor_index].serial.sw->begin(baud_rate);

    return 0;
}

//triggers the sensor and reads the distance in mm, returns -1 on failure, the return result is in result
int32_t read_sensor(int sensor_index) {
    if (sensor_index < 0 || sensor_index >= 4)
        return -1;

    // flush stale data before triggering
    while (sensor_available(sensors[sensor_index].serial))
        sensor_read(sensors[sensor_index].serial);

    sensor_write(sensors[sensor_index].serial, sensors[sensor_index].trigger);

    // wait for 4 bytes: 0xFF + H_byte + L_byte + checksum
    unsigned long start = millis();
    while (sensor_available(sensors[sensor_index].serial) < 4) {
        if ((millis() - start) > MAX_WAIT_TIME_MS)
            return -1;
    }

    uint8_t buffer[4];
    if (sensor_read_bytes(sensors[sensor_index].serial, buffer, 4) < 4)
        return -1;

    if (buffer[0] != 0xFF)
        return -1;

    uint8_t high_byte = buffer[1], low_byte  = buffer[2],
        checksum  = high_byte + low_byte;

    if (checksum != buffer[3])
        return -1;

    float result = (float)((high_byte << 8) | low_byte);
    result = sensors[sensor_index].correct_error(result);
    return ceil(result);
}