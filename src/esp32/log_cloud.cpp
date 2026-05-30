#include "esp32/log_cloud.h"
#include "esp32/config.h"
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <WiFi.h>
#include <PubSubClient.h>

#ifndef WIFI_SSID
#error "Define WIFI_SSID in esp32/config.h"
#endif
#ifndef WIFI_PASSWORD
#error "Define WIFI_PASSWORD in esp32/config.h"
#endif
#ifndef MQTT_BROKER_IP
#error "Define MQTT_BROKER_IP in esp32/config.h"
#endif
#ifndef MQTT_BROKER_PORT
#define MQTT_BROKER_PORT 1883
#endif
#ifndef MQTT_TOPIC
#define MQTT_TOPIC "rtls/position"
#endif
#ifndef MQTT_CLIENT_ID
#define MQTT_CLIENT_ID "rtls-esp32"
#endif

#define WIFI_TIMEOUT_MS     15000
#define MQTT_RETRY_DELAY_MS 2000
#define MQTT_MAX_RETRIES    5

static WiFiClient   wifi_client;
static PubSubClient mqtt_client(wifi_client);

static int wifi_connect() {
    if (WiFi.status() == WL_CONNECTED)
        return 0;

    Serial.printf("[cloud] connecting to WiFi SSID: %s\n", WIFI_SSID);
    WiFi.mode(WIFI_STA);
    WiFi.begin(WIFI_SSID, WIFI_PASSWORD);

    unsigned long start = millis();
    while (WiFi.status() != WL_CONNECTED) {
        if (millis() - start > WIFI_TIMEOUT_MS) {
            Serial.println("[cloud] WiFi connection timed out");
            return -1;
        }
        vTaskDelay(pdMS_TO_TICKS(250));
    }
    Serial.printf("[cloud] WiFi connected, IP: %s\n", WiFi.localIP().toString().c_str());
    return 0;
}

static int mqtt_connect() {
    if (mqtt_client.connected())
        return 0;

    mqtt_client.setServer(MQTT_BROKER_IP, MQTT_BROKER_PORT);
    mqtt_client.setKeepAlive(15);
    mqtt_client.setBufferSize(256);

    for (int attempt = 0; attempt < MQTT_MAX_RETRIES; attempt++) {
        Serial.printf("[cloud] MQTT connect attempt %d/%d to %s:%d\n",
                      attempt + 1, MQTT_MAX_RETRIES,
                      MQTT_BROKER_IP, MQTT_BROKER_PORT);

        if (mqtt_client.connect(MQTT_CLIENT_ID)) {
            Serial.println("[cloud] MQTT connected");
            return 0;
        }
        Serial.printf("[cloud] MQTT connect failed, rc=%d\n", mqtt_client.state());
        vTaskDelay(pdMS_TO_TICKS(MQTT_RETRY_DELAY_MS));
    }
    return -2;
}

int cloud_begin() {
    int rc = wifi_connect();
    if (rc != 0)
        return rc;

    rc = mqtt_connect();
    return rc;
}

int cloud_publish(uint32_t timestamp, const float pos[2]) {
    if (WiFi.status() != WL_CONNECTED) {
        Serial.println("[cloud] WiFi lost, reconnecting...");
        if (wifi_connect() != 0) return -1;
    }
    if (!mqtt_client.connected()) {
        Serial.println("[cloud] MQTT lost, reconnecting...");
        if (mqtt_connect() != 0) return -1;
    }

    char payload[64];
    snprintf(payload, sizeof(payload),
             "{\"ts\":%lu,\"x\":%.2f,\"y\":%.2f}",
             (unsigned long)timestamp, pos[0], pos[1]);

    if (!mqtt_client.publish(MQTT_TOPIC, payload, false)) {
        Serial.printf("[cloud] MQTT publish failed (state=%d)\n", mqtt_client.state());
        return -1;
    }
    return 0;
}

void cloud_loop() {
    if (mqtt_client.connected())
        mqtt_client.loop();
}