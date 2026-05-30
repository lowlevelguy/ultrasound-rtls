#define CONFIG_ESPLOG_SERIAL
#define RX_PIN 19
#define TX_PIN 22

inline HardwareSerial *esplog_serial = &Serial2;