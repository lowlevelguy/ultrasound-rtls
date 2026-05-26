#include "esp32/log_local.h"
#include <LittleFS.h>


static log_queue_entry_t log_queue[2][QUEUE_SIZE];
static uint32_t queue_index = 0;
static uint32_t buffer_index = 0;

int initialize_log_file(fs::FS *fs) {
    if(!fs->exists(LOG_FILE_PATH)) {
        File file = fs->open(LOG_FILE_PATH, FILE_WRITE);
        if(!file) {
            return -1;
        }
        file.println("timestamp,x,y");
        file.close();
    }
    return 0;
}

// saves the contents of the queue to the log file, returns 0 on success and -1 on failure
int save_to_file(fs::FS *fs, uint32_t buffer_index) {

    File file = fs->open(LOG_FILE_PATH, FILE_APPEND);

    if(!file) {
        return -1;
    }

    for(int i = 0; i < queue_index; ++i) {
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
    queue_index = 0;
    return 0;
}

//saves the position to the queue, returns 0 on success and 1 if the queue is full
int on_read(const float* pos, const unsigned long timestamp) {
    log_queue[buffer_index][queue_index].x = pos[0];
    log_queue[buffer_index][queue_index].y = pos[1];
    log_queue[buffer_index][queue_index].timestamp = timestamp;
    queue_index++;

    if(queue_index >= QUEUE_SIZE) {
        buffer_index = 1 - buffer_index;
        xTaskNotify(logger_handle, 1-buffer_index, eSetValueWithOverwrite);
        return 1;
    }

    return 0;
}

void logger_task(void* file_s) {
    uint32_t full_buff_idx;
    xTaskNotifyWait(0, 0, &full_buff_idx, portMAX_DELAY);
    save_to_file((fs::FS*)file_s, full_buff_idx);
}

