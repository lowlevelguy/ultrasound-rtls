#include <Arduino.h>
#include <LittleFS.h>
#include "esp32/config.h"
#include "esp32/log.h"
#include <WiFi.h>
#include <math.h>


void monitor_task(void *param)
{
    while (1)
    {
        vTaskDelay(pdMS_TO_TICKS(5000));

        Serial.printf("\n----------- system status -----------\n");

        // task stack usage
        Serial.printf("[stack] logger:            %u words free\n", uxTaskGetStackHighWaterMark(logger_handle));
        Serial.printf("[stack] uart_receive:       %u words free\n", uxTaskGetStackHighWaterMark(uart_handle));
        Serial.printf("[stack] get_db_last_entry:  %u words free\n", uxTaskGetStackHighWaterMark(get_db_last_entry_handle));
        Serial.printf("[stack] flush_to_db:        %u words free\n", uxTaskGetStackHighWaterMark(flush_to_db_handle));

        // wifi status
        Serial.printf("[wifi]  status: %s\n", WiFi.status() == WL_CONNECTED ? "connected" : "disconnected");
        if (WiFi.status() == WL_CONNECTED) {
            Serial.printf("[wifi]  IP: %s  RSSI: %d dBm\n", WiFi.localIP().toString().c_str(), WiFi.RSSI());
        }

        // mqtt status
        Serial.printf("[mqtt]  status: %s\n", client_connected() ? "connected" : "disconnected");
        Serial.printf("[mqtt]  last_db_entry: %u\n", last_db_entry);

        // file status
        xSemaphoreTake(file_mutex, portMAX_DELAY);
        File f = LittleFS.open(LOG_FILE_PATH, "r");
        if (!f) {
            Serial.println("[file]  could not open log file");
            xSemaphoreGive(file_mutex);
            continue;
        }
        size_t file_size = f.size();
        f.close();
        xSemaphoreGive(file_mutex);
        Serial.printf("[file]  size: %u bytes  (~%u records)\n", file_size, file_size / sizeof(timepos_t));

        // heap
        Serial.printf("[heap]  free: %u bytes  min ever: %u bytes\n", esp_get_free_heap_size(), esp_get_minimum_free_heap_size());

        Serial.printf("-------------------------------------\n\n");
    }
}
void setup()
{
    Serial.begin(115200);
    esplog_serial->begin(9600, SERIAL_8N1, RX_PIN, TX_PIN);
    LittleFS.begin(true);
    file_mutex = xSemaphoreCreateMutex();
    mqtt_mutex = xSemaphoreCreateMutex();
    
    while(cloud_begin() != 0);

    xTaskCreate(
        flush_to_db,
        "flush_to_db_task",
        4096,
        &LittleFS,
        1,
        &flush_to_db_handle
    );
    xTaskCreate(
        get_db_last_entry,
        "get_db_last_entry_task",
        4096,
        nullptr,
        1,
        &get_db_last_entry_handle
    );
    xTaskCreate(
        logger_task,
        "logger_task",
        4096,
        &LittleFS,
        1,
        &logger_handle);
    xTaskCreate(
        uart_receive,
        "uart_receive_task",
        4096,
        nullptr,
        5,
        &uart_handle);
    xTaskCreate(
        monitor_task,
        "monitor",
        4096,
        nullptr,
        1,
        nullptr);
}

void loop() {}