#include <stdlib.h>

#define IND_3D(i, j, k, N) ((i)*(N)*(N) + (j)*(N) + (k))

// u[z][y][x]
void prob_init(double *u, int H, int W, double start_T)
{
  // Set interior points
  for(int i = 1; i < H - 1; i++)
    for(int j = 1; j < W - 1; j++)
      for(int k = 1; k < W - 1; k++)
        u[IND_3D(i, j, k, W)] = start_T;

  // Set boundary conditions
  for(int i = 1; i < H - 1; i++) // Z
    for(int j = 0; j < W - 1; j++) { // X
      u[IND_3D(i, 0    , j, W)] = 0;
      u[IND_3D(i, W - 1, j, W)] = 20;
    }

  for(int i = 1; i < H - 1; i++) // Z
    for(int j = 0; j < W - 1; j++) { // Y
      u[IND_3D(i, j, 0    , W)] = 20;
      u[IND_3D(i, j, W - 1, W)] = 20;
    }

  for(int i = 1; i < W - 1; i++) // Y
    for(int j = 0; j < W - 1; j++) { // X
        u[IND_3D(0    , i, j, W)] = 20;
        u[IND_3D(H - 1, i, j, W)] = 20;
    }

}

void source_init(double *source, int H, int W)
{
  double delta = (2. / W),
    deltasq = delta * delta;
  int
    lx = 0,
    ux = W * 5 / 16,
    ly = 0,
    uy = W / 4,
    lz = H / 6,
    uz = H / 2;

  for(int i = 0; i < H; i++)
    for(int j = 0; j < W; j++)
      for(int k = 0; k < W; k++)
        source[IND_3D(i, j, k, W)] = 0;

  // Place radiator
  for(int i = lz; i < uz; i++)
    for(int j = ly; j < uy; j++)
      for(int k = lx; k < ux; k++)
        source[IND_3D(i, j, k, W)] = 200 * deltasq;

}
