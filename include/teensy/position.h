#pragma once

#include <stdint.h>

/** @brief Performs trilateration to estimate the position of the target.
  *
  * @details
  * Performs trilateration to find the coordinates of the target given its
  * distances, provided in the buffer dist in mm, to three non-collinear anchor
  * points.
  *
  * Algorithm Description:
  * If E denotes the 2x2 matrix of edges, with (k-1)th row given by the vector
  * operation (anchor_pos[0] - anchor_pos[k]), the target point x satisfies the
  * equation Ex = y under ideal circumstances, where y is
  * the 2x1 vector defined as follows:
  *      y[k-1] = (dist[0]**2 - dist[k]**2 + norm_sq(a[k]) - norm_sq(a[0])) / 2
  * We use norm_sq(u) to denote the norm squared of a vector u, and a[k] as a
  * shorthand for anchor_pos[k].
  * This property allows us to recover x when y is known, via the formula:
  *      x = E^{-1} y
  *
  * @param dist: array containing the distances from the target to each anchor
  * in mm.
  * @param pos: array that will contain the estimated position of the target,
  * each coordinate in mm.
  * 
  * @return
  * Returns 0 on success, and -1 if E is non-invertible, meaning that the anchors
  * are collinear.
*/
int position_trilateration(const uint16_t* dist, float* pos);

/** @brief Performs linear ordinary least squares to estimate the position of
  * the target.
  *
  * @details
  * Performs linear ordinary least squares to find the coordinates of the target
  * given its distances, provided in the buffer dist in mm, to four anchor
  * points, three of which are non-collinear.
  * 
  * Algorithm Description:
  * If E denotes the 3x2 matrix of edges, with (k-1)th row given by the vector
  * operation (anchor_pos[0] - anchor_pos[k]), the target point x is the unique
  * minimizer of norm_sq(Ex-y) under ideal circumstances, where y is the 3x1
  * vector defined as follows:
  *      y[k-1] = (dist[0]**2 - dist[k]**2 + norm_sq(a[k]) - norm_sq(a[0])) / 2
  * We use norm_sq(u) to denote the norm squared of a vector u, and a[k] as a
  * shorthand for anchor_pos[k].
  * This property allows us to recover x when y is known, via the formula:
  *      x = (E^T E)^{-1} E^T y
  * 
  * @param dist: array containing the distances from the target to each anchor
  * in mm.
  * @param pos: array that will contain the estimated position of the target,
  * each coordinate in mm.
  * 
  * @return
  * Returns 0 on success, and -1 if (E^T E) is non-invertible, meaning
  * that the condition of three of the anchors being non-collinear is unfulfilled.
  */
int position_ols(const uint16_t* dist, float* pos);

/** @brief Performs linear feasible generalized least squares to estimate the
  * position of the target.
  *
  * @details
  * Performs linear feasile generalized least squares to estimate the position
  * the target.
  * 
  * Algorithm Description:
  * If E denotes the 3x2 matrix of edges, with (k-1)th row given by the vector
  * operation (anchor_pos[0] - anchor_pos[k]), the best unbiased linear estimator
  * hat{x} for the target x, assuming Gaussian i.i.d. noise on the distance
  * measurements, is given by the following expression:
  *      hat{x} = (E^T P^{-1} E)^{-1} E^T P^{-1} y
  * where y is the 3x1 vector defined as follows:
  *      y[k-1] = (dist[0]**2 - dist[k]**2 + norm_sq(a[k]) - norm_sq(a[0])) / 2
  * and P is the covariance matrix of the vector y, given by:
  *     P[i,j] = s**4 / 2  + (s * d[0])**2               (i != j)
  *     P[i,i] = s**4 + s**2 * (d[0]**2 + d[i]**2)
  * We use norm_sq(u) to denote the norm squared of a vector u, a[k] as a
  * shorthand for anchor_pos[k]. Moreover, d[k] denotes the true distance
  * between the target and anchor_pos[k], and s**2 the variance of the noise
  * in the distance measurements.
  * 
  * Given that all d[k]**2 values are unknown, we estimate them with dist[k]**2
  * instead, in order to construct an estimator for P, which will be further
  * used to construct an estimator for the target.
  * Of course, the estimator \hat{P} for P is given by:
  *     \hat{P}[i,j] = s**4 / 2  + (s * dist[0])**2               (i != j)
  *     \hat{P}[i,i] = s**4 + s**2 * (dist[0]**2 + dist[i]**2)
  *
  * @param dist: array containing the distances from the target to each anchor
  * in mm.
  * @param pos: array that will contain the estimated position of the target,
  * each coordinate in mm.
  * 
  * @return
  * Returns 0 on success, and -1 if (E^T P^{-1} E) is non-invertible, meaning
  * that the condition of three of the anchors being non-collinear is unfulfilled.
  */
int position_fgls(const uint16_t* dist, const float sigma_sq, float* pos);

/** @brief Numerically estimates the Maximum Likelihood Estimator for the
  * target's position, assuming homoscedasticity.
  *
  * @details
  * Numerically estimates the Maximum Likelihood Estimator for the target's
  * position, assuming homoscedastic noise, using Guass-Newton iterations to
  * minimize the loss function:
  *     L(x) = norm_sq(d(x) - dist)
  * where we use norm_sq(u) to denote the norm squared of the vector u, and d(x)
  * to denote the array of distances from the point x to the anchors:
  *     d(x)[k] = norm(x - anchor_pos[k])
  *
  * @param dist: array containing the distances from the target to each anchor
  * in mm.
  * @param pos: array that will contain the estimated position of the target,
  * each coordinate in mm.
  *
  * @return
  * Returns 0 on success, and -1 if (J^T J) is non-invertible, where J is the
  * Jacobian of the function (d(x) - dist).
  */
int position_mle(const uint16_t* dist, float* pos);

int position_mle3(const uint16_t* dist, float* pos);
