#pragma once

extern const uint32_t ANCHORS[4][2];

/*
 * Performs trilateration to find the coordinates of the target given its
 * distances to three non-collinear anchor points.
 *
 * Algorithm Description:
 * Let E denote the matrix of edges, with (k-1)th row ANCHORS[0] - ANCHORS[k].
 * The target point x satisfies the equation Ex = y, and hence x = E^{-1}y,
 * where y is the vector defined as follows:
 *      y[k-1] = (dist[0]**2 - dist[k]**2 + norm_sq(a[k]) - norm_sq(a[0])) / 2
 * We use norm_sq(u) to denote the norm squared of a vector u, and a[k] as a
 * shorthand for ANCHORS[k].
 *
 * @param uint64_t* dist: array of distances from the target to the anchors
 * @param int64_t* pos: array that will contain the estimated position of the
 *                      target at the end of function execution.
*/
void position_trilateration(uint64_t* dist, int64_t* pos);

void position_ols(uint64_t* dist, uint64_t* pos);
void position_fgls(uint64_t* dist, uint64_t* pos);