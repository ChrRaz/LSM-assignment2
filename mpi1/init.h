#ifndef INIT_H
#define INIT_H

void prob_init(double *u, int H, int W, double start_T, int rank, int size);
void source_init(double *source, int H, int W, int rank, int size);

#endif /* INIT_H */
