/* jacobi.c - Poisson problem in 3d
 *
 */
#include <math.h>
#ifdef _TOL_PER_IT
#include <stdio.h>
#endif


// size is outer dimension (whole box)
// source is precomputed delta^2*f
int jacobi(double ***u, double ***old_u, double ***source, int size, int iterations, double *threshold)
{
  const double oneoversix = 1./6., term = (*threshold) * (*threshold);
  double d = INFINITY;

  int it = 0;
  while(d > term && it < iterations) {

    d = 0;

    for(int i = 1; i < size - 1; i++)   // z direction
      for(int j = 1; j < size - 1; j++)   // y direction
        for(int k = 1; k < size - 1; k++) { // x direction
          u[i][j][k] = oneoversix * (old_u[i-1][j][k] + old_u[i+1][j][k] + old_u[i][j-1][k] + old_u[i][j+1][k] + old_u[i][j][k-1] + old_u[i][j][k+1] + source[i][j][k]);
          double diff = u[i][j][k] - old_u[i][j][k];
          d += diff * diff;
        }

#ifdef _TOL_PER_IT
      printf("%5d %10lf\n", it, d);
#endif

    it++;

    double ***tmp = u;
    u = old_u;
    old_u = tmp;
  }

  *threshold = sqrt(d);
  return it;
}
