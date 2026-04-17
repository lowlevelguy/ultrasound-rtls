#include "position.h"

#define EPS 1e-9

/********** HELPER FUNCTIONS **********/

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

void matmat_mult(const double* mat1, const uint8_t mat1_rows, const uint8_t mat1_cols,
                const double* mat2, const uint8_t mat2_cols,
                double* output) {
    for (uint8_t i = 0; i < mat1_rows; i++) {
        for (uint8_t j = 0; j < mat2_cols; j++) {
            output[i * mat2_cols + j] = 0;
            for (uint8_t k = 0; k < mat1_cols; k++)
                output[i * mat2_cols + j] +=
                    mat1[i * mat2_cols + k] * mat2[k * mat2_cols + j];
        }
    }
}

int mat2_inv(const double* mat, double* output) {
    const double* det = mat[0]*mat[3] - mat[1]*mat[2];
    if (det < EPS && det > -EPS)
        return -1;

    output[0] = mat[3] / det;
    output[1] = -mat[1] / det;
    output[2] = -mat[2] / det;
    output[3] = mat[0] / det;
    return 0;
}


/********** LIBRARY FUNCTIONS **********/

void position_trilateration(uint64_t* dist, int32_t* pos) {
    // Use the matrix of cofactors to compute the inverse of the edges matrix E
    static const int64_t edges_det =
        (ANCHORS[0][0] - ANCHORS[1][0])*(ANCHORS[0][1] - ANCHORS[2][1])
        - (ANCHORS[0][0] - ANCHORS[2][0])*(ANCHORS[0][1] - ANCHORS[1][1]);
    static const double edges_inv[] = {
        (double)(ANCHORS[0][1] - ANCHORS[2][1]) / det,
        -(double)(ANCHORS[0][1] - ANCHORS[1][1]) / det,
        -(double)(ANCHORS[0][0] - ANCHORS[2][0]) / det,
        (double)(ANCHORS[0][0] - ANCHORS[1][0]) / det
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

int position_ols(uint64_t* dist, uint64_t* pos) {
    static const double edges[] = {
        ANCHORS[0][0] - ANCHORS[1][0],
        ANCHORS[0][1] - ANCHORS[1][1],
        ANCHORS[0][0] - ANCHORS[2][0],
        ANCHORS[0][1] - ANCHORS[2][1],
        ANCHORS[0][0] - ANCHORS[3][0],
        ANCHORS[0][1] - ANCHORS[3][1],
    }, edges_transpose[] = {
        ANCHORS[0][0] - ANCHORS[1][0],
        ANCHORS[0][0] - ANCHORS[2][0],
        ANCHORS[0][0] - ANCHORS[3][0],
        ANCHORS[0][1] - ANCHORS[1][1],
        ANCHORS[0][1] - ANCHORS[2][1],
        ANCHORS[0][1] - ANCHORS[3][1],
    };

    // Compute the necessary matrices
    double edges_product[2*2], edges_product_inv[2*2], full_matrix[2*3];
    matmat_mult(edges_transpose, 2, 3, edges, 2, edges_product); // E^T E
    if (mat2_inv(edges_product, edges_product_inv) == -1) // (E^T E)^{-1}
        return -1;

    // (E^T E)^{-1} E^T
    matmat_mult(edges_product_inv, 2, 2, edges_transpose, 3, full_matrix);


    // Compute the vectors of coefficients
    static const double norms_sq[] = {
        ANCHORS[0][0]*ANCHORS[0][0] + ANCHORS[0][1]*ANCHORS[0][1],
        ANCHORS[1][0]*ANCHORS[1][0] + ANCHORS[1][1]*ANCHORS[1][1],
        ANCHORS[2][0]*ANCHORS[2][0] + ANCHORS[2][1]*ANCHORS[2][1],
        ANCHORS[3][0]*ANCHORS[3][0] + ANCHORS[3][1]*ANCHORS[3][1],
    };
    const double coeffs[] = {
        (dist[0]*dist[0] - dist[1]*dist[1] + norms_sq[1] - norms_sq[0]) / 2,
        (dist[0]*dist[0] - dist[2]*dist[2] + norms_sq[2] - norms_sq[0]) / 2,
        (dist[0]*dist[0] - dist[3]*dist[3] + norms_sq[3] - norms_sq[0]) / 2
    };

    // Estimate the target as (E^T E)^{-1} E^T y
    matvec_mult(full_matrx, 2, 3, coeffs, pos);
}

void position_fgls(uint64_t* dist, uint64_t* pos) {

}