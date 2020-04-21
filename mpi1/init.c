#include <stdlib.h>

#define IND_3D(i, j, k, N) ((i)*(N)*(N) + (j)*(N) + (k))

// u[z][y][x]
void prob_init(double *u, int size, double start_T)
{
  // Set interior points
  for(int i = 1; i < size - 1; i++)
    for(int j = 1; j < size - 1; j++)
      for(int k = 1; k < size - 1; k++)
        u[IND_3D(i, j, k, size)] = start_T;

  // Set boundary conditions
  for(int i = 1; i < size - 1; i++)
    for(int j = 0; j < size - 1; j++) {
      u[IND_3D(i, 0       , j, size)] = 0;
      u[IND_3D(i, size - 1, j, size)] = 20;
    }

  for(int i = 1; i < size - 1; i++)
    for(int j = 0; j < size - 1; j++) {
      u[IND_3D(i, j, 0       , size)] = 20;
      u[IND_3D(i, j, size - 1, size)] = 20;
    }

  for(int i = 1; i < size - 1; i++)
    for(int j = 0; j < size - 1; j++) {
      u[IND_3D(0       , i, j, size)] = 20;
      u[IND_3D(size - 1, i, j, size)] = 20;
    }

}

void source_init(double *source, int size)
{
  double delta = (2. / size),
    deltasq = delta * delta;
  int
    lx = 0,
    ux = size * 5 / 16,
    ly = 0,
    uy = size / 4,
    lz = size / 6,
    uz = size / 2;

  for(int i = 0; i < size; i++)
    for(int j = 0; j < size; j++)
      for(int k = 0; k < size; k++)
        source[IND_3D(i, j, k, size)] = 0;

  // Place radiator
  for(int i = lz; i < uz; i++)
    for(int j = ly; j < uy; j++)
      for(int k = lx; k < ux; k++)
        source[IND_3D(i, j, k, size)] = 200 * deltasq;

}
