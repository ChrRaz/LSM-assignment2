
// u[z][y][x]
void prob_init(double ***u, int size, double start_T)
{
  // Set interior points
#pragma omp parallel for
  for(int i = 1; i < size - 1; i++)
    for(int j = 1; j < size - 1; j++)
      for(int k = 1; k < size - 1; k++)
        u[i][j][k] = start_T;

  // Set boundary conditions
  for(int i = 1; i < size - 1; i++)
    for(int j = 0; j < size - 1; j++) {
      u[i][0]       [j] = 0;
      u[i][size - 1][j] = 20;
    }

  for(int i = 1; i < size - 1; i++)
    for(int j = 0; j < size - 1; j++) {
      u[i][j][0]        = 20;
      u[i][j][size - 1] = 20;
    }

  for(int i = 1; i < size - 1; i++)
    for(int j = 0; j < size - 1; j++) {
      u[0]       [i][j] = 20;
      u[size - 1][i][j] = 20;
    }

}

void source_init(double ***source, int size)
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

#pragma omp parallel for
  for(int i = 0; i < size; i++)
    for(int j = 0; j < size; j++)
      for(int k = 0; k < size; k++)
        source[i][j][k] = 0;

#pragma omp parallel for
  // Place radiator
  for(int i = lz; i < uz; i++)
    for(int j = ly; j < uy; j++)
      for(int k = lx; k < ux; k++)
        source[i][j][k] = 200 * deltasq;

}
