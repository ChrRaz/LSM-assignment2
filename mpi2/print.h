#ifndef _PRINT_H
#define _PRINT_H

void print_binary(char *fname, int H, int W, double *u);

FILE *open_vtk(const char *fname, int H, int W);

void write_vtk(FILE *file, int H, int W, double *u);

#endif
