#include <Arduino.h>

#include "teensy/position.h"
#include "teensy/config.h"

float max_err, min_err, mean, variance;

float rand_stdnormal();

/* Teensy specific functions for cycle count */
inline void enable_cyccounter() {
    ARM_DEMCR |= ARM_DEMCR_TRCENA;
    ARM_DWT_CTRL |= ARM_DWT_CTRL_CYCCNTENA;
}

inline uint32_t get_cyccount() {
    return ARM_DWT_CYCCNT;
}

void setup() {
    enable_cyccounter();
    Serial.begin(115200);

    // Generate inputs
    const float sigma = 5;
    float true_pos[2], pos[2],
        pos_err[500];
        //cycle_counts[301];
    uint32_t start_cycle;
    uint16_t dists[4];
    for (int i = 0; i < 500; i++) {
        true_pos[0] = (float)rand() / (float)RAND_MAX * 4000;
        true_pos[1] = (float)rand() / (float)RAND_MAX * 4000;

        dists[0] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[0][0]) * (true_pos[0] - anchor_pos[0][0])
                 + (true_pos[1] - anchor_pos[0][1]) * (true_pos[1] - anchor_pos[0][1]));
        dists[0] += (uint16_t)(sigma * rand_stdnormal());

        dists[1] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[1][0]) * (true_pos[0] - anchor_pos[1][0])
                 + (true_pos[1] - anchor_pos[1][1]) * (true_pos[1] - anchor_pos[1][1]));
        dists[1] += (uint16_t)(sigma * rand_stdnormal());

        dists[2] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[2][0]) * (true_pos[0] - anchor_pos[2][0])
                 + (true_pos[1] - anchor_pos[2][1]) * (true_pos[1] - anchor_pos[2][1]));
        dists[2] += (uint16_t)(sigma * rand_stdnormal());

        dists[3] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[3][0]) * (true_pos[0] - anchor_pos[3][0])
                 + (true_pos[1] - anchor_pos[3][1]) * (true_pos[1] - anchor_pos[3][1]));
        dists[3] += (uint16_t)(sigma * rand_stdnormal());

//        position_fgls(dists, sigma*sigma, pos);
        position_fgls(dists, sigma*sigma, pos);
        pos_err[i] = sqrtf((pos[0] - true_pos[0]) * (pos[0] - true_pos[0])
            + (pos[1] - true_pos[1]) * (pos[1] - true_pos[1]));
    }

    // Stats for position estimation accuracy
    max_err = pos_err[0];
    min_err = pos_err[0];
    mean = pos_err[0];
    for (int i = 1; i < 500; i++) {
        if (pos_err[i] > max_err)
            max_err = pos_err[i];
        if (pos_err[i] < min_err)
            min_err = pos_err[i];
        mean += pos_err[i];
    }
    mean /= 500;

    variance = 0;
    for (int i = 0; i < 500; i++) {
        variance += (pos_err[i] - mean) * (pos_err[i] - mean);
    }
    variance /= 500;
}

void loop() {
    delay(1000);
    Serial.printf("Min error:      %.2f mm\n", min_err);
    Serial.printf("Max error:      %.2f mm\n", max_err);
    Serial.printf("Mean error:     %.2f mm\n", mean);
    Serial.printf("Std dev error:  %.2f mm\n", sqrtf(variance));
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