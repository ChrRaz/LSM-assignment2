/* main.c - Poisson problem in 3D
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include "print.h"
#include "init.h"

#include "jacobi.h"
#define NMATS 3

#define N_DEFAULT 102

int main(int argc, char *argv[]) {

  int       H             = N_DEFAULT;
  int       W             = N_DEFAULT;
  int       iter_max      = 1000;
  // Unused
  double    tolerance     = 0.1;
  double    start_T       = 0;
  int       output_type   = 0;
  // char     *output_prefix = "poisson_res";
  char     *output_prefix = argv[0];
  char     *output_ext    = "";
  char      output_filename[FILENAME_MAX];
  double   *u             = NULL;

  MPI_Init(&argc, &argv);

  /* get the paramters from the command line */
  if(argc > 2) H           = atoi(argv[1]) + 2; // grid size
  if(argc > 2) W           = atoi(argv[2]) + 2; // grid size
  if(argc > 3) iter_max    = atoi(argv[3]);     // max. no. of iterations
  if(argc > 4) tolerance   = atof(argv[4]);     // tolerance
  if(argc > 5) start_T     = atof(argv[5]);     // start T for all inner grid points
  if(argc > 6) output_type = atoi(argv[6]);     // ouput type

  int world_size, world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_size != 2) {
    fprintf(stderr, "This version can only run on 2 ranks\n");
    return 1;
  }

  // allocate memory
  if ( (u = malloc((H/2 + 2) * W*W * sizeof(double))) == NULL ) {
    perror("array u: allocation failed");
    exit(-1);
  }

  prob_init(u, H/2 + 2, W, start_T, world_rank, world_size);

  // allocate heat source
  double *source;
  if ( (source = malloc((H/2 + 2) * W*W * sizeof(double))) == NULL ) {
    perror("array source: allocation failed");
    exit(-1);
  }

  source_init(source, H/2 + 2, W, world_rank, world_size);

  // allocate old_u
  double *old_u;
  if ( (old_u = malloc((H/2 + 2) * W*W * sizeof(double))) == NULL ) {
    perror("array old_u: allocation failed");
    exit(-1);
  }

  prob_init(old_u, H/2 + 2, W, start_T, world_rank, world_size);

  double t_begin = omp_get_wtime();

  int iters = jacobi(u, old_u, source, H/2 + 2, W, iter_max, world_rank, world_size);
  // int iters = 0;

  if (iters % 2 == 0) {
    double *tmp = old_u;
    old_u = u;
    u = tmp;
  }
  free(old_u);

  double t_end = omp_get_wtime();

  double *result;
  if ( (result = malloc(H*W*W * sizeof(double))) == NULL ) {
    perror("array source: allocation failed");
    exit(-1);
  }

  MPI_Gather(u + W*W, (H/2)*W*W, MPI_DOUBLE, result, (H/2)*W*W, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  if (world_rank == 0) {
    printf("mpi1: %s (%3dx%3d) %d %lf %.1lf %d", argv[0], H - 2, W - 2, iter_max, tolerance, start_T, output_type);
    printf(" : %10.3lf sec, %8lu KB, %d iterations\n", t_end - t_begin, H*W*W * sizeof(double), iters);

    // printf("%s: %d iterations, |x| = %lf\n", argv[0], iters, tolerance);

    // <

    // dump  results if wanted
    switch(output_type) {
    case 0:
      // no output at all
      break;
    case 3:
      output_ext = ".bin";
      sprintf(output_filename, "%s_%dx%d%s", output_prefix, H-2, W-2, output_ext);
      fprintf(stderr, "Write binary dump to %s\n", output_filename);
      print_binary(output_filename, H, W, u);
      break;
    case 4:
      output_ext = ".vtk";
      sprintf(output_filename, "%s_%dx%d%s", output_prefix, H-2, W-2, output_ext);
      fprintf(stderr, "Write VTK file to %s\n", output_filename);
      print_vtk(output_filename, H, W, result);
      break;
    default:
      fprintf(stderr, "Non-supported output type!\n");
      break;
    }
  }

  //output_ext = ".vtk";
  //sprintf(output_filename, "%s_%dx%d_%d%s", output_prefix, H-2, W-2, world_rank, output_ext);
  //fprintf(stderr, "Write VTK file to %s\n", output_filename);
  //print_vtk(output_filename, H/2+2, W, u);

  MPI_Barrier(MPI_COMM_WORLD);

  // de-allocate memory
  free(source);
  free(u);
  free(result);

  MPI_Finalize();

  return(0);
}
