/* jacobi.c - Poisson problem in 3d
 *
 */
#include <math.h>
#ifdef _TOL_PER_IT
#include <stdio.h>
#endif

#define IND_3D(i, j, k, N) ((i)*(N)*(N) + (j)*(N) + (k))


// size is outer dimension (whole box)
// source is precomputed delta^2*f
int jacobi(double *u, double *old_u, double *source, int H, int W, int iterations, double *threshold)
{
  const double oneoversix = 1./6., term = (*threshold) * (*threshold);
  double d = INFINITY;

  int it = 0;
  while(d > term && it < iterations) {

    d = 0;

    for(int i = 1; i < H - 1; i++)   // z direction
      for(int j = 1; j < W - 1; j++)   // y direction
        for(int k = 1; k < W - 1; k++) { // x direction
          u[IND_3D(i, j, k, W)] = oneoversix * (old_u[IND_3D(i-1, j, k, W)] + old_u[IND_3D(i+1, j, k, W)] + old_u[IND_3D(i, j-1, k, W)] + old_u[IND_3D(i, j+1, k, W)] + old_u[IND_3D(i, j, k-1, W)] + old_u[IND_3D(i, j, k+1, W)] + source[IND_3D(i, j, k, W)]);
          double diff = u[IND_3D(i, j, k, W)] - old_u[IND_3D(i, j, k, W)];
          d += diff * diff;
        }

#ifdef _TOL_PER_IT
      printf("%5d %10lf\n", it, d);
#endif

    it++;

    double *tmp = u;
    u = old_u;
    old_u = tmp;
  }

  *threshold = sqrt(d);
  return it;
}
