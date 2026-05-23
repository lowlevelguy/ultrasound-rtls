#include <string.h>
#include <math.h>
#include "position.h"
#include "config.h"

#ifndef CONFIG_ANCHOR_POS
#error "Please define anchor_pos in config.h"
#endif

#define EPS_DET 1e-6f
#define STEP_GN 1e-2f
#define MAX_ITERATIONS_GN 4


/********** HELPER FUNCTIONS **********/

/** @brief Writes the result of the matrix-vector multiplication of mat and vec
  * to output.
  *
  * @details
  * Multiplies the matrix mat of dimension mat_rows * mat_cols by the vector
  * vec of size mat_cols, and writes the resulting vector of size mat_rows to
  * output.
  *
  * @param mat: matrix to left-multiply by, in row-major order.
  * @param mat_rows: number of rows in mat.
  * @param mat_cols: number of columns in mat.
  * @param vec: vector to right-multiply by.
  * @param output: buffer to contain the resultant vector.
  */
void matvec_mult(const float *mat, const uint8_t mat_rows, const uint8_t mat_cols,
                 const float *vec, float *output) {
    for (uint8_t i = 0; i < mat_rows; i++) {
        output[i] = 0;
        for (uint8_t j = 0; j < mat_cols; j++)
            output[i] += mat[i * mat_cols + j] * vec[j];
    }
}

/** @brief Writes the result of the matrix multiplication of mat1 and mat2
  * to output.
  *
  * @details
  * Multiplies the matrix mat1 of dimension mat1_rows * mat1_cols by the matrix
  * mat2 of dimension mat1_cols * mat2_cols, and writes the resulting matrix of
  * dimension mat1_rows * mat2_cols to output.
  *
  * @param mat1: matrix to left-multiply by, in row-major order.
  * @param mat1_rows: number of rows in mat1.
  * @param mat1_cols: number of columns in mat1.
  * @param mat2: matrix to right-multiply by, in row-major order.
  * @param mat2_cols: number of columns in mat2.
  * @param output: buffer to contain the resultant matrix.
  */
void matmat_mult(const float *mat1, const uint8_t mat1_rows, const uint8_t mat1_cols,
                 const float *mat2, const uint8_t mat2_cols,
                 float *output) {
    for (uint8_t i = 0; i < mat1_rows; i++) {
        for (uint8_t j = 0; j < mat2_cols; j++) {
            output[i * mat2_cols + j] = 0;
            for (uint8_t k = 0; k < mat1_cols; k++)
                output[i * mat2_cols + j] +=
                        mat1[i * mat1_cols + k] * mat2[k * mat2_cols + j];
        }
    }
}

/** @brief Computes the inverse of the 2*2 matrix mat when possible, and writes
  * the inverse to output.
  *
  * @param mat: the 2*2 matrix to invert, in row-major order.
  * @param output: buffer to contain the inverse matrix.
  * @return Returns 0 if the matrix is invertible, -1 if not (i.e., |det| < EPS_DET).
  */
int mat2_inv(const float *mat, float *output) {
    const float det = mat[0] * mat[3] - mat[1] * mat[2];
    if (det < EPS_DET && det > -EPS_DET)
        return -1;

    output[0] = mat[3] / det;
    output[1] = -mat[1] / det;
    output[2] = -mat[2] / det;
    output[3] = mat[0] / det;
    return 0;
}

/** @brief Computes the inverse of the 3*3 matrix mat when possible, and writes
  * the inverse to output.
  *
  * @param mat: the 3*3 matrix to invert, in row-major order.
  * @param output: buffer to contain the inverse matrix.
  * @return Returns 0 if the matrix is invertible, -1 if not (i.e., |det| < EPS_DET).
  */
int mat3_inv(const float *mat, float *output) {
    const float det = mat[0] * (mat[4] * mat[8] - mat[5] * mat[7])
        - mat[1] * (mat[3] * mat[8] - mat[5] * mat[6])
        + mat[2] * (mat[3] * mat[7] - mat[4] * mat[6]);
    if (det < EPS_DET && det > -EPS_DET)
        return -1;

    output[0] = (mat[4]*mat[8] - mat[5]*mat[7]) / det;
    output[1] = (mat[2]*mat[7] - mat[8]*mat[1]) / det;
    output[2] = (mat[1]*mat[5] - mat[4]*mat[2]) / det;
    output[3] = (mat[5]*mat[6] - mat[8]*mat[3]) / det;
    output[4] = (mat[0]*mat[8] - mat[6]*mat[2]) / det;
    output[5] = (mat[2]*mat[3] - mat[5]*mat[0]) / det;
    output[6] = (mat[3]*mat[7] - mat[6]*mat[4]) / det;
    output[7] = (mat[1]*mat[6] - mat[7]*mat[0]) / det;
    output[8] = (mat[0]*mat[4] - mat[3]*mat[1]) / det;

    return 0;
}

/** @brief Computes the point-to-anchors distance errors for a given reference
  * point pos and measured distances meas, and writes the anchor-respective\
  * errors to err.
  *
  * @param pos 2-sized array representing the reference point
  * @param meas 4-sized array of the measured target-anchor distances
  * @param err 4-sized buffer that will contain the array of distance errors
  */
void dist_err(const float* pos, const uint16_t* meas, float* err) {
    err[0] = (pos[0] - anchor_pos[0][0])*(pos[0] - anchor_pos[0][0])
         + (pos[1] - anchor_pos[0][1])*(pos[1] - anchor_pos[0][1]);
    err[0] = sqrtf(err[0]);
    err[0] -= meas[0];

    err[1] = (pos[0] - anchor_pos[1][0])*(pos[0] - anchor_pos[1][0])
         + (pos[1] - anchor_pos[1][1])*(pos[1] - anchor_pos[1][1]);
    err[1] = sqrtf(err[1]);
    err[1] -= meas[1];

    err[2] = (pos[0] - anchor_pos[2][0])*(pos[0] - anchor_pos[2][0])
         + (pos[1] - anchor_pos[2][1])*(pos[1] - anchor_pos[2][1]);
    err[2] = sqrtf(err[2]);
    err[2] -= meas[2];

    err[3] = (pos[0] - anchor_pos[3][0])*(pos[0] - anchor_pos[3][0])
         + (pos[1] - anchor_pos[3][1])*(pos[1] - anchor_pos[3][1]);
    err[3] = sqrtf(err[3]);
    err[3] -= meas[3];
}

/********** LIBRARY FUNCTIONS **********/

int position_trilateration(const uint16_t *dist, float *pos) {
    // Use the matrix of cofactors to compute the inverse of the edges matrix E
    static const float edges_det =
            (anchor_pos[0][0] - anchor_pos[1][0]) * (anchor_pos[0][1] - anchor_pos[2][1])
            - (anchor_pos[0][0] - anchor_pos[2][0]) * (anchor_pos[0][1] - anchor_pos[1][1]);

    if (edges_det < EPS_DET && edges_det > -EPS_DET)
        return -1;

    static const float edges_inv[] = {
        (anchor_pos[2][1] - anchor_pos[0][1]) / edges_det,
        -(anchor_pos[1][1] - anchor_pos[0][1]) / edges_det,
        -(anchor_pos[2][0] - anchor_pos[0][0]) / edges_det,
        (anchor_pos[1][0] - anchor_pos[0][0]) / edges_det
    };

    // Compute the vector of coefficients y
    static const float norms_sq[] = {
        anchor_pos[0][0] * anchor_pos[0][0] + anchor_pos[0][1] * anchor_pos[0][1],
        anchor_pos[1][0] * anchor_pos[1][0] + anchor_pos[1][1] * anchor_pos[1][1],
        anchor_pos[2][0] * anchor_pos[2][0] + anchor_pos[2][1] * anchor_pos[2][1]
    };
    const float coeffs[] = {
        (dist[0] * dist[0] - dist[1] * dist[1] + norms_sq[1] - norms_sq[0]) / 2,
        (dist[0] * dist[0] - dist[2] * dist[2] + norms_sq[2] - norms_sq[0]) / 2
    };

    // Estimate the target as E^{1}y
    matvec_mult(edges_inv, 2, 2, coeffs, pos);
    return 0;
}

int position_ols(const uint16_t *dist, float *pos) {
    static const float edges[] = {
        anchor_pos[1][0] - anchor_pos[0][0],
        anchor_pos[1][1] - anchor_pos[0][1],
        anchor_pos[2][0] - anchor_pos[0][0],
        anchor_pos[2][1] - anchor_pos[0][1],
        anchor_pos[3][0] - anchor_pos[0][0],
        anchor_pos[3][1] - anchor_pos[0][1],
    }, edges_transpose[] = {
        anchor_pos[1][0] - anchor_pos[0][0],
        anchor_pos[2][0] - anchor_pos[0][0],
        anchor_pos[3][0] - anchor_pos[0][0],
        anchor_pos[1][1] - anchor_pos[0][1],
        anchor_pos[2][1] - anchor_pos[0][1],
        anchor_pos[3][1] - anchor_pos[0][1],
    };

    // Compute the necessary matrices
    float edges_product[2 * 2], edges_product_inv[2 * 2], full_matrix[2 * 3];
    matmat_mult(edges_transpose, 2, 3, edges, 2, edges_product); // E^T E
    if (mat2_inv(edges_product, edges_product_inv) == -1) // (E^T E)^{-1}
        return -1;

    // (E^T E)^{-1} E^T
    matmat_mult(edges_product_inv, 2, 2, edges_transpose, 3, full_matrix);


    // Compute the vectors of coefficients
    static const float norms_sq[] = {
        anchor_pos[0][0] * anchor_pos[0][0] + anchor_pos[0][1] * anchor_pos[0][1],
        anchor_pos[1][0] * anchor_pos[1][0] + anchor_pos[1][1] * anchor_pos[1][1],
        anchor_pos[2][0] * anchor_pos[2][0] + anchor_pos[2][1] * anchor_pos[2][1],
        anchor_pos[3][0] * anchor_pos[3][0] + anchor_pos[3][1] * anchor_pos[3][1],
    };
    const float coeffs[] = {
        (dist[0] * dist[0] - dist[1] * dist[1] + norms_sq[1] - norms_sq[0]) / 2,
        (dist[0] * dist[0] - dist[2] * dist[2] + norms_sq[2] - norms_sq[0]) / 2,
        (dist[0] * dist[0] - dist[3] * dist[3] + norms_sq[3] - norms_sq[0]) / 2
    };

    // Estimate the target as (E^T E)^{-1} E^T y
    matvec_mult(full_matrix, 2, 3, coeffs, pos);
    return 0;
}

int position_fgls(const uint16_t *dist, const float sigma_sq, float *pos) {
    static const float edges[] = {
        anchor_pos[1][0] - anchor_pos[0][0],
        anchor_pos[1][1] - anchor_pos[0][1],
        anchor_pos[2][0] - anchor_pos[0][0],
        anchor_pos[2][1] - anchor_pos[0][1],
        anchor_pos[3][0] - anchor_pos[0][0],
        anchor_pos[3][1] - anchor_pos[0][1],
    }, edges_transpose[] = {
        anchor_pos[1][0] - anchor_pos[0][0],
        anchor_pos[2][0] - anchor_pos[0][0],
        anchor_pos[3][0] - anchor_pos[0][0],
        anchor_pos[1][1] - anchor_pos[0][1],
        anchor_pos[2][1] - anchor_pos[0][1],
        anchor_pos[3][1] - anchor_pos[0][1],
    };

    if (position_ols(dist, pos) == -1)
        return -1;

    // Estimate P^{-1}
    const float cov[] = {
        sigma_sq*sigma_sq + sigma_sq*(dist[0]*dist[0] + dist[1]*dist[1]),
        sigma_sq*sigma_sq/2 + sigma_sq*dist[0]*dist[0],
        sigma_sq*sigma_sq/2 + sigma_sq*dist[0]*dist[0],

        sigma_sq*sigma_sq/2 + sigma_sq*dist[0]*dist[0],
        sigma_sq*sigma_sq + sigma_sq*(dist[0]*dist[0] + dist[2]*dist[2]),
        sigma_sq*sigma_sq/2 + sigma_sq*dist[0]*dist[0],

        sigma_sq*sigma_sq/2 + sigma_sq*dist[0]*dist[0],
        sigma_sq*sigma_sq/2 + sigma_sq*dist[0]*dist[0],
        sigma_sq*sigma_sq + sigma_sq*(dist[0]*dist[0] + dist[3]*dist[3]),
    };
    float cov_inv[3*3];
    mat3_inv(cov, cov_inv);

    // Compute E^T P^{-1}, (E^T P^{-1} E)^{-1},
    // and the final matrix (E^T P^{-1} E)^{-1} E^T P^{-1}
    float edges_transpose_cov[2*3], weighted_edges_prod[2*2], weighted_edges_prod_inv[2*2], full_matrix[2*3];
    matmat_mult(edges_transpose, 2, 3, cov_inv, 3, edges_transpose_cov);
    matmat_mult(edges_transpose_cov, 2, 3, edges, 2, weighted_edges_prod);
    if (mat2_inv(weighted_edges_prod, weighted_edges_prod_inv) == -1)
            return -1;
    matmat_mult(weighted_edges_prod_inv, 2, 2, edges_transpose_cov, 3, full_matrix);

    // Compute the vectors of coefficients
    static const float norms_sq[] = {
        anchor_pos[0][0] * anchor_pos[0][0] + anchor_pos[0][1] * anchor_pos[0][1],
        anchor_pos[1][0] * anchor_pos[1][0] + anchor_pos[1][1] * anchor_pos[1][1],
        anchor_pos[2][0] * anchor_pos[2][0] + anchor_pos[2][1] * anchor_pos[2][1],
        anchor_pos[3][0] * anchor_pos[3][0] + anchor_pos[3][1] * anchor_pos[3][1],
    };
    const float coeffs[] = {
        (dist[0] * dist[0] - dist[1] * dist[1] + norms_sq[1] - norms_sq[0]) / 2,
        (dist[0] * dist[0] - dist[2] * dist[2] + norms_sq[2] - norms_sq[0]) / 2,
        (dist[0] * dist[0] - dist[3] * dist[3] + norms_sq[3] - norms_sq[0]) / 2
    };

    // Estimate the target as (E^T P E)^{-1} E^T P y
    matvec_mult(full_matrix, 2, 3, coeffs, pos);
    return 0;
}

int position_mle(const uint16_t* dist, float* pos) {
    // Perform Gauss-Newton to find the point with minimal squared error
    float point[2] = {0}, point_shift[2],
        val[4], val_shift[4], val_tf[4],
        J[4*2], Jt[2*4], JtJ[2*2], JtJ_inv[2*2], full_matrix[2*4],
        norm_sq;

    dist_err(point, dist, val);
    for (int i = 0; i < MAX_ITERATIONS_GN; i++) {
        // Compute Jacobian and its transpose
        point_shift[0] = point[0] + STEP_GN;
        point_shift[1] = point[1];
        dist_err(point_shift, dist, val_shift);
        J[0] = (val_shift[0] - val[0]) / STEP_GN;
        J[2] = (val_shift[1] - val[1]) / STEP_GN;
        J[4] = (val_shift[2] - val[2]) / STEP_GN;
        J[6] = (val_shift[3] - val[3]) / STEP_GN;
        Jt[0] = (val_shift[0] - val[0]) / STEP_GN;
        Jt[1] = (val_shift[1] - val[1]) / STEP_GN;
        Jt[2] = (val_shift[2] - val[2]) / STEP_GN;
        Jt[3] = (val_shift[3] - val[3]) / STEP_GN;

        point_shift[0] = point[0];
        point_shift[1] = point[1] + STEP_GN;
        dist_err(point_shift, dist, val_shift);
        J[1] = (val_shift[0] - val[0]) / STEP_GN;
        J[3] = (val_shift[1] - val[1]) / STEP_GN;
        J[5] = (val_shift[2] - val[2]) / STEP_GN;
        J[7] = (val_shift[3] - val[3]) / STEP_GN;
        Jt[4] = (val_shift[0] - val[0]) / STEP_GN;
        Jt[5] = (val_shift[1] - val[1]) / STEP_GN;
        Jt[6] = (val_shift[2] - val[2]) / STEP_GN;
        Jt[7] = (val_shift[3] - val[3]) / STEP_GN;

        // Compute (J^T J)^{-1}
        matmat_mult(Jt, 2, 4, J, 2, JtJ);
        if (mat2_inv(JtJ, JtJ_inv) == -1)
            return -1;

        // Compute (J^T J)^{-1} J^T and right-multiply it by f(point)
        matmat_mult(JtJ_inv, 2, 2, Jt, 4, full_matrix);
        matvec_mult(full_matrix, 2, 4, val, val_tf);

        // Update point
        point[0] -= val_tf[0];
        point[1] -= val_tf[1];

        dist_err(point, dist, val);
    }

    pos[0] = point[0];
    pos[1] = point[1];
    return 0;
}
