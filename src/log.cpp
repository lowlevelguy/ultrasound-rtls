#include "log.h"
#include <LittleFS.h>

static log_queue_entry_t log_queue[QUEUE_SIZE];
static int queue_index = 0;

int initialize_log_file(fs::FS &fs) {
    if(!fs.exists(LOG_FILE_PATH)) {
        File file = fs.open(LOG_FILE_PATH, FILE_WRITE);
        if(!file) {
            return -1;
        }
        file.println("timestamp,x,y");
        file.close();
    }
    return 0;
}

// saves the contents of the queue to the log file, returns 0 on success and -1 on failure
int save_to_file(fs::FS &fs) {

    File file = fs.open(LOG_FILE_PATH, FILE_APPEND);

    if(!file) {
        return -1;
    }

    for(int i = 0; i < queue_index; ++i) {
        if(!file.printf("%lu,%.2f,%.2f\n",log_queue[i].timestamp , log_queue[i].x, log_queue[i].y)) {
            file.close();
            return -1;
        }
    }
    file.close();
    queue_index = 0;
    return 0;
}

//saves the position to the queue, returns 0 on success and 1 if the queue is full
int save_to_queue(const float* pos, const unsigned long timestamp) {
    if(queue_index >= QUEUE_SIZE) {
        return 1;
    }

    log_queue[queue_index].x = pos[0];
    log_queue[queue_index].y = pos[1];
    log_queue[queue_index].timestamp = timestamp;
    queue_index++;

    return 0;
}