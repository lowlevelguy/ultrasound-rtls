#include "position.h"
#include <math.h>

#define EPS 1e-6f

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

/** @brief Computes the inverse of the 3*3 matrix mat when possible, and writes
  * the inverse to output.
  *
  * @param mat: the 3*3 matrix to invert, in row-major order.
  * @param output: buffer to contain the inverse matrix.
  * @return Returns 0 if the matrix is invertible, -1 if not (i.e., |det| < EPS).
  */
int mat2_inv(const float *mat, float *output) {
    const float det = mat[0] * mat[3] - mat[1] * mat[2];
    if (det < EPS && det > -EPS)
        return -1;

    output[0] = mat[3] / det;
    output[1] = -mat[1] / det;
    output[2] = -mat[2] / det;
    output[3] = mat[0] / det;
    return 0;
}

int mat3_inv(const float *mat, float *output) {
    const float det = mat[0] * (mat[4] * mat[8] - mat[5] * mat[7])
        - mat[1] * (mat[3] * mat[8] - mat[5] * mat[6])
        + mat[2] * (mat[3] * mat[7] - mat[4] * mat[6]);
    if (det < EPS && det > -EPS)
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


/********** LIBRARY FUNCTIONS **********/

int position_trilateration(const uint16_t *dist, float *pos) {
    // Use the matrix of cofactors to compute the inverse of the edges matrix E
    static const float edges_det =
            (anchor_pos[0][0] - anchor_pos[1][0]) * (anchor_pos[0][1] - anchor_pos[2][1])
            - (anchor_pos[0][0] - anchor_pos[2][0]) * (anchor_pos[0][1] - anchor_pos[1][1]);

    if (edges_det < EPS && edges_det > -EPS)
        return -1;

    static const float edges_inv[] = {
        (anchor_pos[0][1] - anchor_pos[2][1]) / edges_det,
        -(anchor_pos[0][1] - anchor_pos[1][1]) / edges_det,
        -(anchor_pos[0][0] - anchor_pos[2][0]) / edges_det,
        (anchor_pos[0][0] - anchor_pos[1][0]) / edges_det
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
        anchor_pos[0][0] - anchor_pos[1][0],
        anchor_pos[0][1] - anchor_pos[1][1],
        anchor_pos[0][0] - anchor_pos[2][0],
        anchor_pos[0][1] - anchor_pos[2][1],
        anchor_pos[0][0] - anchor_pos[3][0],
        anchor_pos[0][1] - anchor_pos[3][1],
    }, edges_transpose[] = {
        anchor_pos[0][0] - anchor_pos[1][0],
        anchor_pos[0][0] - anchor_pos[2][0],
        anchor_pos[0][0] - anchor_pos[3][0],
        anchor_pos[0][1] - anchor_pos[1][1],
        anchor_pos[0][1] - anchor_pos[2][1],
        anchor_pos[0][1] - anchor_pos[3][1],
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

int position_fgls(const uint16_t *dist, float *pos) {
    static const float edges[] = {
        anchor_pos[0][0] - anchor_pos[1][0],
        anchor_pos[0][1] - anchor_pos[1][1],
        anchor_pos[0][0] - anchor_pos[2][0],
        anchor_pos[0][1] - anchor_pos[2][1],
        anchor_pos[0][0] - anchor_pos[3][0],
        anchor_pos[0][1] - anchor_pos[3][1],
    }, edges_transpose[] = {
        anchor_pos[0][0] - anchor_pos[1][0],
        anchor_pos[0][0] - anchor_pos[2][0],
        anchor_pos[0][0] - anchor_pos[3][0],
        anchor_pos[0][1] - anchor_pos[1][1],
        anchor_pos[0][1] - anchor_pos[2][1],
        anchor_pos[0][1] - anchor_pos[3][1],
    };

    if (position_ols(dist, pos) == -1)
        return -1;

    // Compute hat{s}**2
    float sigma_sq = 0;
    for (int i = 0; i < 3; i++)
        sigma_sq += dist[i] - sqrtf(
            (pos[0] - anchor_pos[i][0])*(pos[0] - anchor_pos[i][0])
            + (pos[1] - anchor_pos[i][1])*(pos[1] - anchor_pos[i][1]));
    sigma_sq /= 4;

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
    matmat_mult(edges_transpose, 2, 3, cov, 3, edges_transpose_cov);
    matmat_mult(edges_transpose_cov, 2, 3, edges, 2, weighted_edges_prod);
    if (mat2_inv(weighted_edges_prod, weighted_edges_prod_inv) == -1)
            return -1;
    matmat_mult(weighted_edges_prod_inv, 2, 2, edges_transpose_cov, 2, full_matrix);

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
