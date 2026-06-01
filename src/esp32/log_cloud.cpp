#include "esp32/log.h"
#include "esp32/config.h"
#include <freertos/FreeRTOS.h>
#include <freertos/task.h>
#include <LittleFS.h>
#include "FS.h"
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

#define WIFI_TIMEOUT_MS     10000
#define MQTT_RETRY_DELAY_MS 2000
#define MQTT_MAX_RETRIES 5
#define MAX_TP_SIZE 20
#define LOG_FILE_PATH "/history.log"


uint32_t last_db_entry = 0;
static volatile bool sent_request = false;
static volatile bool received_response = false;
static WiFiClient   wifi_client;
static PubSubClient mqtt_client(wifi_client);
TaskHandle_t flush_to_db_handle = NULL;
TaskHandle_t get_db_last_entry_handle = NULL;
SemaphoreHandle_t mqtt_mutex = NULL;

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
    mqtt_client.setCallback(callback);

    for (int attempt = 0; attempt < MQTT_MAX_RETRIES; attempt++) {
        Serial.printf("[cloud] MQTT connect attempt %d/%d to %s:%d\n",
                      attempt + 1, MQTT_MAX_RETRIES,
                      MQTT_BROKER_IP, MQTT_BROKER_PORT);

        if (mqtt_client.connect(MQTT_CLIENT_ID) && mqtt_client.subscribe(MQTT_TOPIC_RESPONSE)) {
            Serial.println("[cloud] MQTT connected and subscribed");
            return 0;
        }
        Serial.printf("[cloud] MQTT connect or subscribe failed, rc=%d\n", mqtt_client.state());
        vTaskDelay(pdMS_TO_TICKS(MQTT_RETRY_DELAY_MS));
    }
    return -1;
}

int cloud_begin() {
    int rc = wifi_connect();
    if (rc != 0)
        return rc;

    rc = mqtt_connect();
    return rc;
}

int  cloud_publish(timepos_t* tp, size_t tp_size) {
    if(tp_size > MAX_TP_SIZE) return -1;
    int rc = cloud_begin();
    if (rc != 0) return rc;


    if (!mqtt_client.publish(MQTT_TOPIC, (uint8_t*)tp, tp_size*sizeof(timepos_t))) {
        Serial.printf("[cloud] MQTT publish failed (state=%d)\n", mqtt_client.state());
        return -1;
    }
    return 0;
}

bool client_connected() {
    return mqtt_client.connected();
}

void cloud_loop() {
    if (mqtt_client.connected())
        mqtt_client.loop();
}

void callback(char* topic, byte* message, unsigned int length) {
    if(!sent_request) return;

    if(strcmp(topic, "rtls/position/response") == 0){
        if(length != sizeof(uint32_t)) {
            return;
        }
        memcpy(&last_db_entry, message, sizeof(uint32_t));
        received_response = true;
    }
}
void flush_to_db(void* param) {
    uint32_t time_stamp;
    timepos_t buffer[MAX_TP_SIZE];
    fs::FS* file_s = (fs::FS*)param;
    size_t curr_pos;
    while(1) {
        xTaskNotifyWait(0, 0, &time_stamp, portMAX_DELAY);
        xSemaphoreTake(file_mutex, portMAX_DELAY);
        File file = file_s->open(LOG_FILE_PATH, FILE_READ);
        if(file.size() == 0) {
            file.close();
            xSemaphoreGive(file_mutex);
            continue;
        }
        curr_pos = file.size();
        size_t to_flush = 0;
        uint32_t curr_timestamp;
        
        while(curr_pos>sizeof(timepos_t)) {
            curr_pos -= sizeof(timepos_t);
            file.seek(curr_pos, SeekSet);
            file.readBytes((char*)&curr_timestamp, sizeof(uint32_t));
            to_flush += 1;
            if(curr_timestamp == last_db_entry) {
                break;
            }
        }
        file.close();
        xSemaphoreGive(file_mutex);
        
        size_t loops_num = ((to_flush) / MAX_TP_SIZE);
        size_t rest_to_flush = to_flush % MAX_TP_SIZE;
        
        for(int i = 0; i < loops_num; ++i) {

            xSemaphoreTake(file_mutex, portMAX_DELAY);
            file = file_s->open(LOG_FILE_PATH, FILE_READ);
            file.seek(curr_pos);
            file.readBytes((char*)buffer, MAX_TP_SIZE * sizeof(timepos_t));
            file.close();
            xSemaphoreGive(file_mutex);
            curr_pos += MAX_TP_SIZE * sizeof(timepos_t);
            xSemaphoreTake(mqtt_mutex, portMAX_DELAY);
            cloud_publish(buffer, MAX_TP_SIZE);
            xSemaphoreGive(mqtt_mutex);
        }
        if(rest_to_flush > 0) {
            xSemaphoreTake(file_mutex, portMAX_DELAY);
            file = file_s->open(LOG_FILE_PATH, FILE_READ);
            file.seek(curr_pos);
            file.readBytes((char*)buffer, rest_to_flush*sizeof(timepos_t));
            file.close();
            xSemaphoreGive(file_mutex);
        }
        xSemaphoreTake(mqtt_mutex, portMAX_DELAY);
        cloud_publish(buffer, rest_to_flush);
        xSemaphoreGive(mqtt_mutex);
    }
}

void get_db_last_entry(void* param) {
    uint8_t request[1] = {1};
    while(1) {
        vTaskDelay(pdMS_TO_TICKS(DB_FLUSH_PEORIOD));

        xSemaphoreTake(mqtt_mutex, portMAX_DELAY);
        bool published = mqtt_client.publish(MQTT_TOPIC_REQUEST, request, false);
        xSemaphoreGive(mqtt_mutex);
        if(!published) {
            Serial.printf("failed to publish request to the database\n");
            continue;
        }
        
        sent_request = true;
        uint32_t start_time = millis();
        while(!received_response) {
            if(millis() - start_time > MAX_RESPONSE_TIMEOUT) {
                break;
            }
            xSemaphoreTake(mqtt_mutex, portMAX_DELAY);
            cloud_loop();
            xSemaphoreGive(mqtt_mutex);
            vTaskDelay(1);
        }
        sent_request = false;  
        received_response = false;
        xTaskNotify(flush_to_db_handle, last_db_entry, eSetValueWithOverwrite);
    }
}