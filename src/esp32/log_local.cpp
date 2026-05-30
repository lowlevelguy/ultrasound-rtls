#include "esp32/log_local.h"
#include "esp32/log_cloud.h"
#include "esp32/config.h"
#include <atomic>

static log_queue_entry_t log_queue[2][QUEUE_SIZE];
static uint32_t queue_index = 0;
static std::atomic<uint32_t> buffer_index = 0;
TaskHandle_t logger_handle = NULL;
TaskHandle_t uart_handle = NULL;
static uint32_t bad_pkts_counter = 0;
static QueueHandle_t _cloud_queue = NULL;

void set_cloud_queue(QueueHandle_t q) {
    _cloud_queue = q;
}

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
static int save_to_file(fs::FS *file_s, uint32_t buff_idx) {

    File file = file_s->open(LOG_FILE_PATH, FILE_APPEND);

    if(!file) {
        return -1;
    }

    for(int i = 0; i < QUEUE_SIZE; ++i) {
        if(!file.printf("%lu,%.2f,%.2f\n",
            log_queue[buff_idx][i].timestamp, 
            log_queue[buff_idx][i].x, 
            log_queue[buff_idx][i].y)) 
        {
            file.close();
            return -1;
        }
    }
    file.close();
    return 0;
}

static uint8_t checksum8(uint8_t* buf, uint8_t buf_sz) {
    uint8_t cs = 0;
    for (uint8_t i = 0; i < buf_sz; i++)
        cs += buf[i];

    return cs;
}

//saves the position to the queue, sends notification to logger task if queue is full
static void on_read(const float* pos, const uint32_t timestamp) {
    uint32_t buff_idx = buffer_index.load();
    log_queue[buff_idx][queue_index].x = pos[0];
    log_queue[buff_idx][queue_index].y = pos[1];
    log_queue[buff_idx][queue_index].timestamp = timestamp;
    queue_index++;

    if(queue_index >= QUEUE_SIZE) {
        uint32_t full_buff_index = buffer_index.load();
        buffer_index.fetch_xor(1);
        xTaskNotify(logger_handle, full_buff_index, eSetValueWithOverwrite);
        queue_index = 0;
    }
}

void log_local_push(const float* pos, uint32_t timestamp) {
    on_read(pos, timestamp);
}

void uart_receive(void* param) {
    uint8_t buffer[sizeof(log_packet_t)];
    log_packet_t log_packet;
    ack_packet_t ack;
    ack.start = START_BYTE;
    while(1) {
        while(esplog_serial->available() < (int)sizeof(log_packet_t)) {
            vTaskDelay(pdMS_TO_TICKS(100));
        }
        
        esplog_serial->readBytes(buffer, sizeof(log_packet_t));
        memcpy(&log_packet, buffer, sizeof(log_packet_t));

        if((log_packet.start != START_BYTE) || 
            (checksum8((uint8_t*)&log_packet, sizeof(log_packet_t) - sizeof(uint8_t)) != log_packet.checksum)) {
            Serial.printf("bad packet recieved number: %u\n", bad_pkts_counter++);
            continue;
        }
        
        on_read(log_packet.pos, log_packet.timestamp);

        if (_cloud_queue != NULL) {
            log_queue_entry_t entry = {
                .x         = log_packet.pos[0],
                .y         = log_packet.pos[1],
                .timestamp = (unsigned long)log_packet.timestamp
            };
            xQueueSend(_cloud_queue, &entry, 0);
        }

        ack.timestamp = log_packet.timestamp;
        esplog_serial->write((uint8_t*)&ack, sizeof(ack_packet_t));
    }
}

void logger_task(void* file_s) {
    uint32_t full_buff_idx;

    while(1){
        xTaskNotifyWait(0, 0, &full_buff_idx, portMAX_DELAY);
        save_to_file((fs::FS*)file_s, full_buff_idx);
    }
}