/* jacobi.c - Poisson problem in 3d
 *
 */
#include <stdio.h>
#include <mpi.h>

#include <math.h>

#define IND_3D(i, j, k, W) ((i)*(W)*(W) + (j)*(W) + (k))


// size is outer dimension (whole box)
// source is precomputed delta^2*f
int jacobi(double *u, double *old_u, double *source, int H, int W, int iterations, int rank, int size)
{

  const double oneoversix = 1./6.;

  int it = 0;
  while(it < iterations) {

    if (rank == 0) {
      MPI_Send(&old_u[IND_3D(H-2, 0, 0, W)], W*W, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
      MPI_Recv(&old_u[IND_3D(H-1, 0, 0, W)], W*W, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (rank > 0 && rank < size - 1) {
      MPI_Recv(&old_u[IND_3D(  0, 0, 0, W)], W*W, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&old_u[IND_3D(  1, 0, 0, W)], W*W, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
      MPI_Send(&old_u[IND_3D(H-2, 0, 0, W)], W*W, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
      MPI_Recv(&old_u[IND_3D(H-1, 0, 0, W)], W*W, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (rank == size - 1) {
      MPI_Recv(&old_u[IND_3D(  0, 0, 0, W)], W*W, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(&old_u[IND_3D(  1, 0, 0, W)], W*W, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
    }

    for(int i = 1; i < H - 1; i++)   // z direction
      for(int j = 1; j < W - 1; j++)   // y direction
        for(int k = 1; k < W - 1; k++) { // x direction
          u[IND_3D(i, j, k, W)] = oneoversix * (old_u[IND_3D(i-1, j, k, W)] + old_u[IND_3D(i+1, j, k, W)] + old_u[IND_3D(i, j-1, k, W)] + old_u[IND_3D(i, j+1, k, W)] + old_u[IND_3D(i, j, k-1, W)] + old_u[IND_3D(i, j, k+1, W)] + source[IND_3D(i, j, k, W)]);
        }

    it++;

    double *tmp = u;
    u = old_u;
    old_u = tmp;
  }

  return it;
}
