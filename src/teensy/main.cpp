#include <Arduino.h>

#include "teensy/config.h"
#include "teensy/log.h"
#include "teensy/position.h"

#define SENSOR_READ_DELAY_MS 10

uint16_t dists[4];
float pos[2];
uint32_t timestamp = 0;
void setup() {
    Serial.begin(115200);

    sensor_begin(0, 9600);
    sensor_begin(1, 9600);
    sensor_begin(2, 9600);
    sensor_begin(3, 9600);

    tlog_serial->begin(9600);
}

void loop() {
    int32_t temp;
    if ((temp = read_sensor(0)) < 0) {
        Serial.println("Error reading sensor 0.");
        return;
    }
    dists[0] = temp;
    delay(SENSOR_READ_DELAY_MS);

    if ((temp = read_sensor(1)) < 0) {
        Serial.println("Error reading sensor 1.");
        return;
    }
    dists[1] = temp;
    delay(SENSOR_READ_DELAY_MS);

    if ((temp = read_sensor(2)) < 0) {
        Serial.println("Error reading sensor 2.");
        return;
    }
    dists[2] = temp;
    delay(SENSOR_READ_DELAY_MS);

    if ((temp = read_sensor(3)) < 0) {
        Serial.println("Error reading sensor 3.");
        return;
    }
    dists[3] = temp;
    delay(SENSOR_READ_DELAY_MS);

    if (position_mle(dists, pos) < 0) {
        Serial.println("Error computing position with MLE");
        return;
    }

    Serial.printf("Distances: %u, %u, %u, %u.\n", dists[0], dists[1], dists[2], dists[3]);
    Serial.printf("Estimated position: %f, %f. Error from true position: %f.\n",
        pos[0], pos[1], sqrtf((pos[0] - 750) * (pos[0] - 750) + (pos[1] - 750) * (pos[1] - 750)));
    log_timepos(timestamp, pos);

    timestamp++;
    delay(500);
}