#ifndef READER_H
#define READER_H

#include <petscmat.h>

extern Mat read_petsc_matrix(const char *filename, MPI_Comm comm);

#endif 

