#include "sensor.h"

int init_sensor(int sensor_index, int uart_index, int baud_rate, int RX, int TX) {
    //input validation
    if(sensor_index < 0 || sensor_index >= 4) {
        return 0;
    }
    SENSORS[sensor_index] = new HardwareSerial(uart_index);
    // initialize the sensor at the given index 
    //with the specified baud rate and RX/TX pins
    if (SENSORS[sensor_index] == NULL) {
        return 0;
    }

    SENSORS[sensor_index]->begin(baud_rate, SERIAL_8N1, RX, TX);
    return 1;
}


//triggers the sensor and reads the distance
uint16_t read_sensor(int sensor_index) {
    //input validation 
    if (sensor_index < 0 || sensor_index >= 4) {
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
    SENSORS[sensor_index]->write(SENSOR_TRIGGER_BYTE);

    //sensor has multiple modes of operation, we chose to use the mode 4
    // in this mode , when triggered the sensor sends 32 bits of data
    // in the format 0xff + H_byte + L_byte + checksum
    // where H_byte and L_byte are the high and low bytes of the sensor reading
    // and checksum is the lower 8 bits of the sum of H_byte and L_byte
    unsigned long start = millis();
    while (SENSORS[sensor_index]->available() < 4) {
        if((millis() - start) > MAX_WAIT_TIME) {
            return -1;
        }
    }
    
    SENSORS[sensor_index]->readBytes(buffer, 4);

    if(buffer[0] == 0xff) {
        uint8_t high_bytes = buffer[1];
        uint8_t low_bytes = buffer[2];

        uint8_t checksum = (uint8_t)(((uint16_t)low_bytes + high_bytes) & 0xFF);
        
        if(checksum != buffer[3]) {
            return -1;
        }
        uint16_t distance = (uint16_t)(high_bytes << 8) | low_bytes;
        return distance; 
    
    }
    return -1;
}