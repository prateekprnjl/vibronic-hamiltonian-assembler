#include <petscmat.h>
#include "reader.h"

Mat read_petsc_matrix(const char *filename, MPI_Comm comm)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    Mat A;

    ierr = PetscViewerBinaryOpen(comm, filename, FILE_MODE_READ, &viewer);
    if (ierr) {
        PetscPrintf(comm, "Error: Could not open PETSc file %s\n", filename);
        return NULL;
    }

    ierr = MatCreate(comm, &A); CHKERRABORT(comm, ierr);
    ierr = MatSetFromOptions(A); CHKERRABORT(comm, ierr);
    ierr = MatLoad(A, viewer); CHKERRABORT(comm, ierr);

    PetscViewerDestroy(&viewer);

    return A;
}

