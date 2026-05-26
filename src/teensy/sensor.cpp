#include "teensy/sensor.h"

#define MAX_WAIT_TIME_MS 100
#define SENSOR_TRIGGER_BYTE 0x01

static int sensor_available(SerialHandle_t s) {
    return (s.type == SerialHandle_t::HW)
        ? s.hw->available()
        : s.sw->available();
}

static int sensor_read(SerialHandle_t s) {
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

static void sensor_read_bytes(SerialHandle_t s, uint8_t* buffer, size_t len) {
    if (s.type == SerialHandle_t::HW)
        s.hw->readBytes(buffer, len);
    else
        s.sw->readBytes(buffer, len);
}

int init_sensor(int sensor_index, const SerialHandle_t* handle, int baud_rate, int RX, int TX) {
    if (sensor_index < 0 || sensor_index >= 4)
        return -1;
    if (handle == NULL)
        return -1;

    memcpy(&SENSORS[sensor_index], handle, sizeof(SerialHandle_t));

    if (handle->type == SerialHandle_t::HW)
        handle->hw->begin(baud_rate, SERIAL_8N1, RX, TX);
    else
        handle->sw->begin(baud_rate); // RX/TX already set in constructor in main

    return 0;
}

//triggers the sensor and reads the distance in mm, returns -1 on failure, the return result is in result
int32_t read_sensor(int sensor_index) {
    
    if (sensor_index < 0 || sensor_index >= 4)
        return -1;

    // flush stale data before triggering
    while (sensor_available(SENSORS[sensor_index]))
        sensor_read(SENSORS[sensor_index]);

    sensor_write(SENSORS[sensor_index], SENSOR_TRIGGER_BYTE);

    // wait for 4 bytes: 0xFF + H_byte + L_byte + checksum
    unsigned long start = millis();
    while (sensor_available(SENSORS[sensor_index]) < 4) {
        if ((millis() - start) > MAX_WAIT_TIME_MS)
            return -1;
    }

    uint8_t buffer[4];
    sensor_read_bytes(SENSORS[sensor_index], buffer, 4);

    if (buffer[0] != 0xFF)
        return -1;

    uint8_t high_byte = buffer[1], low_byte  = buffer[2],
        checksum  = (uint8_t)(((uint16_t)high_byte + low_byte) & 0xFF);

    if (checksum != buffer[3])
        return -1;

    float result = (float)((high_byte << 8) | low_byte);
    result = (result + 21.5) / 0.9896; // correcting the measurement error using linear interpolation
    return (int32_t)result;
}