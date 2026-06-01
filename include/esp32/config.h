#pragma once

#define CONFIG_ESPLOG_SERIAL
#define RX_PIN 19
#define TX_PIN 22

inline HardwareSerial *esplog_serial = &Serial2;

#define WIFI_SSID        "TUNISIETELECOM-2.4G-M4c6"
#define WIFI_PASSWORD    "Ec6YgkGm"
#define MAX_RESPONSE_TIMEOUT 2000
#define MQTT_BROKER_IP   "192.168.100.14"
#define MQTT_BROKER_PORT 1883
#define MQTT_TOPIC       "rtls/position/raw"
#define MQTT_TOPIC_REQUEST "rtls/position/request"
#define MQTT_TOPIC_RESPONSE "rtls/position/response"
#define DB_FLUSH_PEORIOD 30000
#define MQTT_CLIENT_ID   "rtls-esp32"