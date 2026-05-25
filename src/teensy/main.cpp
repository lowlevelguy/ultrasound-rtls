#include <Arduino.h>

#include "teensy/position.h"
#include "teensy/config.h"

float max_count, min_count, mean, variance;

float rand_stdnormal();

void setup() {
    Serial.begin(115200);

    // Generate inputs
    const float sigma = 5;
    float true_pos[2], pos[301][2],
        //pos_err[500];
        cycle_counts[301];
    uint32_t start_cycle;
    uint16_t dists[301][4];
    for (int i = 0; i < 301; i++) {
        true_pos[0] = (float)rand() / (float)RAND_MAX * 4000;
        true_pos[1] = (float)rand() / (float)RAND_MAX * 4000;

        dists[i][0] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[0][0]) * (true_pos[0] - anchor_pos[0][0])
                 + (true_pos[1] - anchor_pos[0][1]) * (true_pos[1] - anchor_pos[0][1]));
        dists[i][0] += (uint16_t)(sigma * rand_stdnormal());

        dists[i][1] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[1][0]) * (true_pos[0] - anchor_pos[1][0])
                 + (true_pos[1] - anchor_pos[1][1]) * (true_pos[1] - anchor_pos[1][1]));
        dists[i][1] += (uint16_t)(sigma * rand_stdnormal());

        dists[i][2] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[2][0]) * (true_pos[0] - anchor_pos[2][0])
                 + (true_pos[1] - anchor_pos[2][1]) * (true_pos[1] - anchor_pos[2][1]));
        dists[i][2] += (uint16_t)(sigma * rand_stdnormal());

        dists[i][3] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[3][0]) * (true_pos[0] - anchor_pos[3][0])
                 + (true_pos[1] - anchor_pos[3][1]) * (true_pos[1] - anchor_pos[3][1]));
        dists[i][3] += (uint16_t)(sigma * rand_stdnormal());

//        position_fgls(dists, sigma*sigma, pos);
//        position_mle(dists, pos)
    }

    for (int i = 0; i < 301; i++) {
        portDISABLE_INTERRUPTS();
        start_cycle = esp_cpu_get_ccount();
        position_mle(dists[i], pos[i]);
        cycle_counts[i] = (float)(esp_cpu_get_ccount() - start_cycle);
        portENABLE_INTERRUPTS();
    }

    // Stats for position estimation accuracy
    max_count = cycle_counts[1];
    min_count = cycle_counts[1];
    mean = cycle_counts[1];
    for (int i = 2; i < 301; i++) {
        if (cycle_counts[i] > max_count)
            max_count = cycle_counts[i];
        if (cycle_counts[i] < min_count)
            min_count = cycle_counts[i];
        mean += cycle_counts[i];
    }
    mean /= 300;

    variance = 0;
    for (int i = 1; i < 301; i++) {
        variance += (cycle_counts[i] - mean) * (cycle_counts[i] - mean);
    }
    variance /= 300;
}

void loop() {
    delay(1000);
    Serial.printf("Min cycle count:      %.2f mm\n", min_count);
    Serial.printf("Max cycle count:      %.2f mm\n", max_count);
    Serial.printf("Mean cycle count:     %.2f mm\n", mean);
    Serial.printf("Std dev cycle count:  %.2f mm\n", sqrtf(variance));
}

float rand_stdnormal() {
    static int hasSpare = 0;
    static float spare;

    if (hasSpare) {
        hasSpare = 0;
        return spare;
    }

    hasSpare = 1;

    float u, v, s;
    do {
        u = ((float)rand() / ((float) RAND_MAX)) * 2.0 - 1.0;
        v = ((float)rand() / ((float) RAND_MAX)) * 2.0 - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    s = sqrtf(-2.0 * logf(s) / s);
    spare = v * s;

    return u * s;
}