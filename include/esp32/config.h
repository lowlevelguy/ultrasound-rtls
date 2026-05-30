#pragma once

#define CONFIG_ESPLOG_SERIAL
#define RX_PIN 19
#define TX_PIN 22

inline HardwareSerial *esplog_serial = &Serial2;

#define WIFI_SSID        "wifi_name"
#define WIFI_PASSWORD    "password"
#define MQTT_BROKER_IP   "192.168.1.131"
#define MQTT_BROKER_PORT 1883
#define MQTT_TOPIC       "rtls/position"
#define MQTT_CLIENT_ID   "rtls-esp32"