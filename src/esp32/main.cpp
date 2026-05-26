#include <Arduino.h>
#include <LittleFS.h>

#include "esp32/log_local.h"

unsigned long last_write = 0;
int counter = 0;

void print_log_file() {

    File file = LittleFS.open(LOG_FILE_PATH, FILE_READ);

    if(!file) {
        Serial.println("Failed to open log file for reading");
        return;
    }

    Serial.println("\n===== LOG FILE CONTENT =====");

    while(file.available()) {
        Serial.write(file.read());
    }

    Serial.println("===== END OF FILE =====\n");

    file.close();
}

void setup() {

    Serial.begin(115200);
    delay(2000);

    Serial.println("Starting LittleFS test...");

    if(!LittleFS.begin(true)) {
        Serial.println("LittleFS mount failed");
        return;
    }

    Serial.println("LittleFS mounted");

    if(initialize_log_file(LittleFS) != 0) {
        Serial.println("Failed to initialize log file");
        return;
    }

    Serial.println("Log file initialized");
}

void loop() {

    if(millis() - last_write >= 1000) {

        last_write = millis();

        float pos[2];

        pos[0] = counter * 1.1f;
        pos[1] = counter * 2.2f;

        if(save_to_queue(pos, millis()) == 0) {

            Serial.printf(
                "Queued: x=%.2f y=%.2f t=%lu\n",
                pos[0],
                pos[1],
                millis()
            );

        } else {

            Serial.println("Queue full");
        }

        counter++;

        // save every 5 entries
        if(counter % 5 == 0) {

            if(save_to_file(LittleFS) == 0) {

                Serial.println("Saved queue to file");

                print_log_file();

            } else {

                Serial.println("Failed to save queue");
            }
        }
    }
}