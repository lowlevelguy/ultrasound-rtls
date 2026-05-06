#include <Arduino.h>

#include "position.h"

const float anchor_pos[4][2] = {
    {0, 0},
    {400, 0},
    {0, 400},
    {400, 400}
};

float max_pos_err, min_pos_err, mean_pos_err, variance_pos_err;

float rand_stdnormal();

void setup() {
    Serial.begin(115200);

    // Generate inputs
    const float sigma = 5;
    float true_pos[2], pos[2], pos_err[500];
    uint16_t dists[4];
    for (int i = 0; i < 500; i++) {
        true_pos[0] = (float)rand() / (float)RAND_MAX * 400;
        true_pos[1] = (float)rand() / (float)RAND_MAX * 400;

        dists[0] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[0][0]) * (true_pos[0] - anchor_pos[0][0])
                 + (true_pos[1] - anchor_pos[0][1]) * (true_pos[1] - anchor_pos[0][1]));
//        dists[0] += (uint16_t)(sigma * rand_stdnormal());

        dists[1] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[1][0]) * (true_pos[0] - anchor_pos[1][0])
                 + (true_pos[1] - anchor_pos[1][1]) * (true_pos[1] - anchor_pos[1][1]));
//        dists[1] += (uint16_t)(sigma * rand_stdnormal());

        dists[2] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[2][0]) * (true_pos[0] - anchor_pos[2][0])
                 + (true_pos[1] - anchor_pos[2][1]) * (true_pos[1] - anchor_pos[2][1]));
//        dists[2] += (uint16_t)(sigma * rand_stdnormal());

        dists[3] = (uint16_t)sqrtf((true_pos[0] - anchor_pos[3][0]) * (true_pos[0] - anchor_pos[3][0])
                 + (true_pos[1] - anchor_pos[3][1]) * (true_pos[1] - anchor_pos[3][1]));
//        dists[3] += (uint16_t)(sigma * rand_stdnormal());

        position_fgls(dists, pos);
        pos_err[i] = sqrtf((pos[0] - true_pos[0]) * (pos[0] - true_pos[0])
                        + (pos[1] - true_pos[1]) * (pos[1] - true_pos[1]));
    }

    // Stats for position estimation accuracy
    max_pos_err = pos_err[0];
    min_pos_err = pos_err[0];
    mean_pos_err = pos_err[0];
    for (int i = 1; i < 500; i++) {
        if (pos_err[i] > max_pos_err)
            max_pos_err = pos_err[i];
        if (pos_err[i] < min_pos_err)
            min_pos_err = pos_err[i];
        mean_pos_err += pos_err[i];
    }
    mean_pos_err /= 500;

    variance_pos_err = 0;
    for (int i = 0; i < 500; i++) {
        variance_pos_err += (pos_err[i] - mean_pos_err) * (pos_err[i] - mean_pos_err);
    }
    variance_pos_err /= 500;
}

void loop() {
    delay(1000);
    Serial.printf("Min position error:      %.2f mm\n", min_pos_err);
    Serial.printf("Max position error:      %.2f mm\n", max_pos_err);
    Serial.printf("Mean position error:     %.2f mm\n", mean_pos_err);
    Serial.printf("Std dev position error:  %.2f mm\n", sqrtf(variance_pos_err));
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