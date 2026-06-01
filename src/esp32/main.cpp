#include <Arduino.h>
#include <LittleFS.h>
#include "esp32/config.h"
#include "esp32/log.h"
#include <PubSubClient.h>
#include <WiFi.h>
#include<math.h>
bool received_response = false;
bool sent_request = false;
WiFiClient   wifi_client;
PubSubClient mqtt_client(wifi_client);

void monitor_task(void* param) {
    while(1) {
        vTaskDelay(pdMS_TO_TICKS(5000));

        Serial.printf("----------- stack high water marks: ----------\n");
        Serial.printf("logger task: %u words free\n", uxTaskGetStackHighWaterMark(logger_handle));
        Serial.printf("uart_receive task: %u words free\n", uxTaskGetStackHighWaterMark(uart_handle));
        File f = LittleFS.open(LOG_FILE_PATH, "r");
        if (!f) {
            Serial.println("[monitor] could not open log file");
            continue;
        }
        f.close();
    }
}

void setup() {
    Serial.begin(115200);
    esplog_serial->begin(9600, SERIAL_8N1, RX_PIN, TX_PIN);
    LittleFS.begin(true);

    xTaskCreate(
        logger_task,
        "logger_task",
        4096,
        &LittleFS,
        1,
        &logger_handle
    );
    xTaskCreate(
        uart_receive,
        "uart_receive_task",
        4096,
        nullptr,
        5,
        &uart_handle
    );
    xTaskCreate(
        monitor_task,
        "monitor", 
        4096, 
        nullptr, 
        1, 
        nullptr
    );
}
