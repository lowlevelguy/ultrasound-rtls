#include <Arduino.h>

#include "teensy/log.h"
#include "teensy/config.h"

#ifndef CONFIG_TLOG_SERIAL
#error "Please define tlog_serial in config.h"
#endif

#define WRITE_TIMEOUT_MS 10
#define ACK_TIMEOUT_MS 100

uint8_t checksum8(uint8_t* buf, uint8_t buf_sz) {
    uint8_t cs = 0;
    for (uint8_t i = 0; i < buf_sz; i++)
        cs += buf[i];

    return cs;
}

int log_timepos(uint32_t timestamp, float pos[2]) {
    log_packet_t log_packet;
    ack_packet_t ack_packet;
    log_packet.start = START_BYTE;
    log_packet.timestamp = timestamp;
    log_packet.pos[0] = pos[0];
    log_packet.pos[1] = pos[1];
    log_packet.checksum =
        checksum8((uint8_t*)&log_packet, sizeof(log_packet_t) - sizeof(uint8_t));

    // Wait for the interface to be available for writing
    unsigned long start = millis();
    while (tlog_serial->availableForWrite() < sizeof(log_packet_t)) {
        if (millis() - start > WRITE_TIMEOUT_MS)
            return ERROR_LOG_WRITE_TIMEOUT;
    }

    if (tlog_serial->write((uint8_t*)&log_packet, sizeof(log_packet_t))
        != sizeof(log_packet_t))
        return ERROR_LOG_WRITE_FAILURE;

    start = millis();
    while (tlog_serial->available() < sizeof(ack_packet_t)) {
        if (millis() - start > ACK_TIMEOUT_MS)
            return ERROR_LOG_ACK_TIMEOUT;
    }

    if (tlog_serial->readBytes((uint8_t*)&ack_packet, sizeof(ack_packet_t))
        != sizeof(ack_packet_t))
        return ERROR_LOG_READ_FAILURE;

    if (ack_packet.start != START_BYTE || ack_packet.timestamp != timestamp) {
        Serial.printf("Ack: %x %u\n", ack_packet.start, ack_packet.timestamp);
        return ERROR_LOG_ACK_INVALID;
    }
    return 0;
}
