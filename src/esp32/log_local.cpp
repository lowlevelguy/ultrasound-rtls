#include "esp32/log.h"
#include "esp32/config.h"
#include <atomic>

static timepos_t log_queue[2][QUEUE_SIZE];
static std::atomic<uint32_t> queue_index = 0;
static std::atomic<uint32_t> curr_queue = 0;
TaskHandle_t logger_handle = NULL;
TaskHandle_t uart_handle = NULL;
SemaphoreHandle_t file_mutex = NULL;
static uint32_t bad_pkts_counter = 0;


// saves the contents of the queue to the log file, returns 0 on success and -1 on failure
static int save_to_file(fs::FS *file_s, uint32_t curr_full_queue) {

    File file = file_s->open(LOG_FILE_PATH, FILE_APPEND);

    if(!file) {
        Serial.printf("could not open the file");
        return -1;
    }

    if(!file.write((uint8_t*)log_queue[curr_full_queue], sizeof(log_queue[0]))) {
        file.close();
        return -1;
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
static void log_local_push(const float* pos, const uint32_t timestamp) {
    uint32_t curr_queue_tmp = curr_queue.load();
    uint32_t idx = queue_index.fetch_add(1);
    log_queue[curr_queue_tmp][idx].pos[0] = pos[0];
    log_queue[curr_queue_tmp][idx].pos[1] = pos[1];
    log_queue[curr_queue_tmp][idx].timestamp = timestamp;

    Serial.printf("updated queue\n");
    if(idx+1 >= QUEUE_SIZE) {
        Serial.printf("queue full\n");
        uint32_t curr_full_queue = curr_queue.load();
        queue_index.store(0);
        curr_queue.fetch_xor(1);
        xTaskNotify(logger_handle, curr_full_queue, eSetValueWithOverwrite);
    }
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
        Serial.print("received packet\n");
        
        esplog_serial->readBytes(buffer, sizeof(log_packet_t));
        memcpy(&log_packet, buffer, sizeof(log_packet_t));

        if((log_packet.start != START_BYTE) || 
            (checksum8((uint8_t*)&log_packet, sizeof(log_packet_t) - sizeof(uint8_t)) != log_packet.checksum)) {
            Serial.printf("bad packet recieved number: %u\n", bad_pkts_counter++);
            continue;
        }
        Serial.printf("position: %f, %f\n", log_packet.pos[0], log_packet.pos[1]);
        log_local_push(log_packet.pos, log_packet.timestamp);

        ack.timestamp = log_packet.timestamp;
        esplog_serial->write((uint8_t*)&ack, sizeof(ack_packet_t));
        Serial.printf("sent ack: %x %u\n", ack.start, ack.timestamp);
    }
}

void logger_task(void* file_s) {
    uint32_t curr_full_queue;

    while(1){
        xTaskNotifyWait(0, 0, &curr_full_queue, portMAX_DELAY);
        //lock
        xSemaphoreTake(file_mutex, portMAX_DELAY);
        save_to_file((fs::FS*)file_s, curr_full_queue);
        //unlock
        xSemaphoreGive(file_mutex);
    }
}