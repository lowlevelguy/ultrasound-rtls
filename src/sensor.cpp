#include "sensor.h"
#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
HardwareSerial *SENSORS[MAX_SENSORS];


int init_sensor(int sensor_index, int uart_index, int baud_rate, int RX, int TX) {
    //input validation
    if(sensor_index < 0 || sensor_index >= MAX_SENSORS) {
        return 0;
    }

    if(SENSORS[sensor_index] != NULL) {
        delete SENSORS[sensor_index];
        SENSORS[sensor_index] = NULL;
    }

    SENSORS[sensor_index] = new HardwareSerial(uart_index);
    // initialize the sensor at the given index 
    //with the specified baud rate and RX/TX pins

    SENSORS[sensor_index]->begin(baud_rate, SERIAL_8N1, RX, TX);
    return 1;
}


//triggers the sensor and reads the distance
int32_t read_sensor(int sensor_index) {
    //input validation 
    if (sensor_index < 0 || sensor_index >= MAX_SENSORS) {
        return -1; 
    }

    //initialize buffer for sensor reading
    uint8_t buffer[4];
    if(SENSORS[sensor_index] == NULL) {
        return -1;
    }

    //trigger the sensor
    //according to the sheet of the sensor model AJ-SR04M, you have to send
    //the value 0x01 to its uart port to trigger it to send its measurement
    while (SENSORS[sensor_index]->available()) {
        SENSORS[sensor_index]->read();
    }
    SENSORS[sensor_index]->write(SENSOR_TRIGGER_BYTE);

    //sensor has multiple modes of operation, we chose to use the mode 4
    // in this mode , when triggered the sensor sends 32 bits of data
    // in the format 0xff + H_byte + L_byte + checksum
    // where H_byte and L_byte are the high and low bytes of the sensor reading
    // and checksum is the lower 8 bits of the sum of H_byte and L_byte
    unsigned long start = millis();
    while (SENSORS[sensor_index]->available() < 4) {
        if((millis() - start) > MAX_WAIT_TIME_MS) {
            return -1;
        }
        vTaskDelay(1);
    }
    
    SENSORS[sensor_index]->readBytes(buffer, 4);

    if(buffer[0] == 0xFF) {
        uint8_t high_byte = buffer[1];
        uint8_t low_byte = buffer[2];

        uint8_t checksum = (uint8_t)(((uint16_t)low_byte + high_byte) & 0xFF);
        
        if(checksum != buffer[3]) {
            return -1;
        }
        int32_t distance = (int32_t)((high_byte << 8) | low_byte);
        return distance; 
    
    }
    return -1;
}