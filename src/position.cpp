#include "position.h"

/* Multiplies a matrix of dimension mat_rows * mat_cols by a vector of size
 * mat_cols, and writes the resulting vector of size mat_rows to the output
 * buffer.
 */
void matvec_mult(const double* mat, const uint8_t mat_rows, const uint8_t mat_cols,
                const double* vec, double* output) {
    for (uint8_t i = 0; i < mat_rows; i++) {
        output[i] = 0;
        for (uint8_t j = 0; j < mat_cols; j++)
            output[i] += mat[i * mat_cols + j] * vec[j];
    }
}

void position_trilateration(uint64_t* dist, int32_t* pos) {
    // Use the matrix of cofactors to compute the inverse of the edges matrix E
    static const int64_t edges_det =
        (ANCHORS[0][0] - ANCHORS[1][0])*(ANCHORS[0][1] - ANCHORS[2][1])
        - (ANCHORS[0][0] - ANCHORS[2][0])*(ANCHORS[0][1] - ANCHORS[1][1]);
    static const double edges_inv[] = {
        (1.0/det) * (ANCHORS[0][1] - ANCHORS[2][1]),
        (-1.0/det) * (ANCHORS[0][1] - ANCHORS[1][1]),
        (-1.0/det) * (ANCHORS[0][0] - ANCHORS[2][0]),
        (1.0/det) * (ANCHORS[0][0] - ANCHORS[1][0])
    };

    // Compute the vector of coefficients y
    static const double norms_sq[] = {
        ANCHORS[0][0]*ANCHORS[0][0] + ANCHORS[0][1]*ANCHORS[0][1],
        ANCHORS[1][0]*ANCHORS[1][0] + ANCHORS[1][1]*ANCHORS[1][1],
        ANCHORS[2][0]*ANCHORS[2][0] + ANCHORS[2][1]*ANCHORS[2][1]
    };
    const double coeffs[] = {
        (dist[0]*dist[0] - dist[1]*dist[1] + norms_sq[1] - norms_sq[0]) / 2,
        (dist[0]*dist[0] - dist[2]*dist[2] + norms_sq[2] - norms_sq[0]) / 2
    };

    // Estimate the target as E^{1}y
    matvec_mult(edges_inv, 2, 2, coeffs, pos);
}

void position_ols(uint64_t* dist, uint64_t* pos) {

}

void position_fgls(uint64_t* dist, uint64_t* pos) {

}