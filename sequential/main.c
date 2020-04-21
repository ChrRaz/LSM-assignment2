/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "print.h"
#include "init.h"

#include "jacobi.h"
#define NMATS 3

#define N_DEFAULT 102

int main(int argc, char *argv[]) {

  int       N             = N_DEFAULT;
  int       iter_max      = 1000;
  double    tolerance     = 0.1;
  double    start_T       = 0;
  int       output_type   = 0;
  // char     *output_prefix = "poisson_res";
  char     *output_prefix = argv[0];
  char     *output_ext    = "";
  char      output_filename[FILENAME_MAX];
  double   *u             = NULL;


  /* get the paramters from the command line */
  if(argc > 1) N           = atoi(argv[1]) + 2; // grid size
  if(argc > 2) iter_max    = atoi(argv[2]);     // max. no. of iterations
  if(argc > 3) tolerance   = atof(argv[3]);     // tolerance
  if(argc > 4) start_T     = atof(argv[4]);     // start T for all inner grid points
  if(argc > 5) output_type = atoi(argv[5]);     // ouput type

  double outtol = tolerance;

  // allocate memory
  if ( (u = malloc(N * N * N * sizeof(double))) == NULL ) {
    perror("array u: allocation failed");
    exit(-1);
  }

  prob_init(u, N, start_T);

  // allocate heat source
  double *source;
  if ( (source = malloc(N * N * N * sizeof(double))) == NULL ) {
    perror("array source: allocation failed");
    exit(-1);
  }

  source_init(source, N);

  // allocate old_u
  double *old_u;
  if ( (old_u = malloc(N * N * N * sizeof(double))) == NULL ) {
    perror("array old_u: allocation failed");
    exit(-1);
  }

  prob_init(old_u, N, start_T);

  double t_begin = omp_get_wtime();

  int iters = jacobi(u, old_u, source, N, iter_max, &outtol);

  if (iters % 2 == 0) {
    double *tmp = old_u;
    old_u = u;
    u = tmp;
  }
  free(old_u);

  double t_end = omp_get_wtime();

  printf("RESULT: %s %3d %d %lf %.1lf %d", argv[0], N - 2, iter_max, tolerance, start_T, output_type);
  printf(" : %10.3lf sec, %8lu KB, %d iterations, |x| = %lf\n", t_end - t_begin, (N * sizeof(double **) + N * N * sizeof(double *) + N * N * N * sizeof(double)) * NMATS / 1024, iters, outtol);

  // printf("%s: %d iterations, |x| = %lf\n", argv[0], iters, tolerance);

  // <

  // dump  results if wanted
  switch(output_type) {
  case 0:
    // no output at all
    break;
  case 3:
    output_ext = ".bin";
    sprintf(output_filename, "%s_%d%s", output_prefix, N-2, output_ext);
    fprintf(stderr, "Write binary dump to %s\n", output_filename);
    print_binary(output_filename, N, u);
    break;
  case 4:
    output_ext = ".vtk";
    sprintf(output_filename, "%s_%d%s", output_prefix, N-2, output_ext);
    fprintf(stderr, "Write VTK file to %s\n", output_filename);
    print_vtk(output_filename, N, u);
    break;
  default:
    fprintf(stderr, "Non-supported output type!\n");
    break;
  }

  // de-allocate memory
  free(source);
  free(u);

  return(0);
}
