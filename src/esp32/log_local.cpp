#include "esp32/log_local.h"
#include <LittleFS.h>
#include "esp32/config.h"

static log_queue_entry_t log_queue[2][QUEUE_SIZE];
static uint32_t queue_index = 0;
static uint32_t buffer_index = 0;

int initialize_log_file(fs::FS *file_s) {
    if(!file_s->exists(LOG_FILE_PATH)) {
        File file = file_s->open(LOG_FILE_PATH, FILE_WRITE);
        if(!file) {
            return -1;
        }
        file.println("timestamp,x,y");
        file.close();
    }
    return 0;
}

// saves the contents of the queue to the log file, returns 0 on success and -1 on failure
int save_to_file(fs::FS *file_s, uint32_t buffer_index) {

    File file = file_s->open(LOG_FILE_PATH, FILE_APPEND);

    if(!file) {
        return -1;
    }

    for(int i = 0; i < QUEUE_SIZE; ++i) {
        if(!file.printf("%lu,%.2f,%.2f\n",
            log_queue[buffer_index][i].timestamp, 
            log_queue[buffer_index][i].x, 
            log_queue[buffer_index][i].y)) 
        {
            file.close();
            return -1;
        }
    }
    file.close();
    return 0;
}

uint8_t checksum8(uint8_t* buf, uint8_t buf_sz) {
    uint8_t cs = 0;
    for (uint8_t i = 0; i < buf_sz; i++)
        cs += buf[i];

    return cs;
}

void uart_recieve(void* param) {
    while(1) {
        while(esplog_serial->available() < sizeof(log_packet_t)) {
            vTaskDelay(pdMS_TO_TICKS(100));
        }
        uint8_t buffer[sizeof(log_packet_t)];
        esplog_serial->readBytes(buffer, sizeof(log_packet_t));
        log_packet_t log_packet;
        memcpy(&log_packet, buffer, sizeof(log_packet_t));

        int* res = (int*) param;
        if((log_packet.start != START_BYTE) || 
            (checksum8((uint8_t*)&log_packet, sizeof(log_packet_t) - sizeof(uint8_t)) != log_packet.checksum)) {
            *res = -1;
            continue;
        }
        
        on_read(log_packet.pos, log_packet.timestamp);

        ack_packet_t ack;
        ack.start = START_BYTE;
        ack.timestamp = log_packet.timestamp;
        esplog_serial->write((uint8_t*)&ack, sizeof(ack_packet_t));
    }
}

//saves the position to the queue, sends notification to logger task if queue is full
void on_read(const float* pos, const uint32_t timestamp) {
    log_queue[buffer_index][queue_index].x = pos[0];
    log_queue[buffer_index][queue_index].y = pos[1];
    log_queue[buffer_index][queue_index].timestamp = timestamp;
    queue_index++;

    if(queue_index >= QUEUE_SIZE) {
        buffer_index = 1 - buffer_index;
        xTaskNotify(logger_handle, 1-buffer_index, eSetValueWithOverwrite);
        queue_index = 0;
    }
}

void logger_task(void* file_s) {
    uint32_t full_buff_idx;

    while(1){
        xTaskNotifyWait(0, 0, &full_buff_idx, portMAX_DELAY);
        save_to_file((fs::FS*)file_s, full_buff_idx);
    }
}