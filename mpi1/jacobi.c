/* jacobi.c - Poisson problem in 3d
 *
 */
#include <mpi.h>

#include <math.h>
#ifdef _TOL_PER_IT
#include <stdio.h>
#endif

#define IND_3D(i, j, k, N) ((i)*(N)*(N) + (j)*(N) + (k))


// size is outer dimension (whole box)
// source is precomputed delta^2*f
int jacobi(double *u, double *old_u, double *source, int N, int iterations, double *threshold, int rank)
{
  const double oneoversix = 1./6., term = (*threshold) * (*threshold);
  double d = INFINITY;

  int it = 0;
  while(d > term && it < iterations) {

    d = 0;

    int startz = rank == 0 ? 1   : N/2;
    int endz   = rank == 0 ? N/2 : N-1;

    for(int i = startz; i < endz; i++)   // z direction
      for(int j = 1; j < N - 1; j++)   // y direction
        for(int k = 1; k < N - 1; k++) { // x direction
          u[IND_3D(i, j, k, N)] = oneoversix * (old_u[IND_3D(i-1, j, k, N)] + old_u[IND_3D(i+1, j, k, N)] + old_u[IND_3D(i, j-1, k, N)] + old_u[IND_3D(i, j+1, k, N)] + old_u[IND_3D(i, j, k-1, N)] + old_u[IND_3D(i, j, k+1, N)] + source[IND_3D(i, j, k, N)]);
          double diff = u[IND_3D(i, j, k, N)] - old_u[IND_3D(i, j, k, N)];
          d += diff * diff;
        }

#ifdef _TOL_PER_IT
      printf("%5d %10lf\n", it, d);
#endif

    if (rank == 0) {
      MPI_Send(u + N*N*N/2 - N*N, N*N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
      MPI_Recv(u + N*N*N/2, N*N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(u + N*N*N/2 - N*N, N*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(u + N*N*N/2, N*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    it++;

    double *tmp = u;
    u = old_u;
    old_u = tmp;
  }

  *threshold = sqrt(d);
  return it;
}
