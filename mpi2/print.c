/* $Id: print.c,v 1.1 2019/12/12 15:03:38 gbarbd Exp gbarbd $ */
#include <stdio.h>
#include <inttypes.h>
#include <mpi.h>
#include "print.h"

#define IND_3D(i, j, k, N) ((i)*(N)*(N) + (j)*(N) + (k))

static int
is_little_endian(void) {
    int num = 1;
    return (*((char *)&num) == 1);
}

void
print_binary(char *fname, int H, int W, double *u) {

    FILE *f_ptr;
    size_t written;
    size_t items = H*W*W;

    if ( (f_ptr = fopen(fname, "w")) == NULL ) {
       perror("No output! fopen()");
       return;
    }

    written = fwrite(u, sizeof(double), items, f_ptr);

    if ( written != items ) {
	    fprintf(stderr, "Writing failed:  only %lu of %lu items saved!\n",
		written, items);
    }

    fclose(f_ptr);
}

FILE *open_vtk(const char *fname, int H, int W) {

    FILE *f_ptr;
    size_t items = H*W*W;
    int b;
    unsigned char tmp;

    if ( (f_ptr = fopen(fname, "w")) == NULL ) {
       perror("No output! fopen()");
       return NULL;
    }

    // Write VTK file header
    fprintf(f_ptr, "# vtk DataFile Version 3.0\n");
    fprintf(f_ptr, "saved from function print_vtk.\n");
    fprintf(f_ptr, "BINARY\n");
    fprintf(f_ptr, "DATASET STRUCTURED_POINTS\n");
    fprintf(f_ptr, "DIMENSIONS %d %d %d\n", W, W, H);
    fprintf(f_ptr, "ORIGIN %d %d %d\n", 0, 0, 0);
    fprintf(f_ptr, "SPACING %d %d %d\n", 1, 1, 1);
    fprintf(f_ptr, "POINT_DATA %lu\n", items);
    fprintf(f_ptr, "SCALARS %s %s 1\n", "gray", "double");
    fprintf(f_ptr, "LOOKUP_TABLE default\n");

    return f_ptr;
}

void write_vtk(FILE *file, int H, int W, double *u)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    size_t written;
    size_t items = H*W*W;
    if ( is_little_endian() ) {
        // System is little endian, so we need to reverse the byte order.
        written = 0;
        for (int i = 0; i < items; ++i) {
            uint64_t crnt = *(uint64_t *)(u + i); // Get double as int

            // Reverse byte order and write to file
            crnt = (crnt & 0x00000000FFFFFFFF) << 32 | (crnt & 0xFFFFFFFF00000000) >> 32;
            crnt = (crnt & 0x0000FFFF0000FFFF) << 16 | (crnt & 0xFFFF0000FFFF0000) >> 16;
            crnt = (crnt & 0x00FF00FF00FF00FF) << 8  | (crnt & 0xFF00FF00FF00FF00) >> 8;
            written += fwrite(&crnt, sizeof(uint64_t), 1, file);
        }
    } else {
        // System is big endian, so just dump the data.
        written = fwrite(u, sizeof(double), items, file);
    }

    if ( written != items ) {
	    fprintf(stderr, "Writing failed:  only %lu of %lu items saved!\n",
		written, items);
    } else {
        fprintf(stderr, "[%d] Wrote %d items\n", rank, written);
    }

}
