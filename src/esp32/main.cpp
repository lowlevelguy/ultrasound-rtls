#include <Arduino.h>
#include <LittleFS.h>
#include "esp32/log_local.h"
#include "esp32/log_cloud.h"
#include "esp32/config.h"

static QueueHandle_t cloud_queue = NULL;

static void cloud_task(void* param) {
    if (cloud_begin() != 0) {
        Serial.println("[cloud] init failed, task stopped");
        vTaskDelete(NULL);
        return;
    }
    log_queue_entry_t entry;
    while (true) {
        if (xQueueReceive(cloud_queue, &entry, pdMS_TO_TICKS(100)) == pdTRUE) {
            const float pos[2] = { entry.x, entry.y };
            if (cloud_publish((uint32_t)entry.timestamp, pos) != 0)
                Serial.println("[cloud] publish failed");
        }
        cloud_loop();
    }
}

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
    initialize_log_file(&LittleFS);

    cloud_queue = xQueueCreate(20, sizeof(log_queue_entry_t));
    set_cloud_queue(cloud_queue);

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
    xTaskCreatePinnedToCore(
        cloud_task,
        "cloud",
        8192,
        nullptr,
        1,
        nullptr,
        0
    );
}


void loop() {}