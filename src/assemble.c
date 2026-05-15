#include <petscmat.h>
#include <mpi.h>
#include <sys/stat.h>

static PetscErrorCode check_file_exists(const char *fname, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
        struct stat buffer;
        PetscCheck(stat(fname, &buffer) == 0, comm, PETSC_ERR_FILE_OPEN,
                   "ERROR: Required block file not found: %s", fname);
    }

    MPI_Barrier(comm);
    return 0;
}

static PetscErrorCode load_block(const char *fname, MPI_Comm comm, Mat *A)
{
    PetscErrorCode ierr;
    PetscViewer viewer;
    int rank;

    MPI_Comm_rank(comm, &rank);

    ierr = check_file_exists(fname, comm); CHKERRQ(ierr);

    if (rank == 0)
        PetscPrintf(comm, "Loading block: %s\n", fname);

    ierr = PetscViewerBinaryOpen(comm, fname, FILE_MODE_READ, &viewer); CHKERRQ(ierr);
    ierr = MatCreate(comm, A); CHKERRQ(ierr);
    ierr = MatSetFromOptions(*A); CHKERRQ(ierr);
    ierr = MatLoad(*A, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    return 0;
}

PetscErrorCode assemble_full_matrix(MPI_Comm comm)
{
    PetscErrorCode ierr;
    int rank;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0)
        PetscPrintf(comm, "\n=== Assembling Full Nested Matrix (4x4 Blocks) ===\n");

    Mat neu_neu, neu_X, neu_A, neu_B;
    Mat X_neu,  X_X,  X_A,  X_B;
    Mat A_neu,  A_X,  A_A,  A_B;
    Mat B_neu,  B_X,  B_A,  B_B;

    load_block("blocks/neu_result.bin",  comm, &neu_neu);
    load_block("blocks/neu_X.bin",       comm, &neu_X);
    load_block("blocks/neu_A.bin",       comm, &neu_A);
    load_block("blocks/neu_B.bin",       comm, &neu_B);

    load_block("blocks/X_neu_conj.bin",  comm, &X_neu);
    load_block("blocks/X_result.bin",    comm, &X_X);
    load_block("blocks/X_A.bin",         comm, &X_A);
    load_block("blocks/X_B.bin",         comm, &X_B);

    load_block("blocks/A_neu_conj.bin",  comm, &A_neu);
    load_block("blocks/A_X_conj.bin",    comm, &A_X);
    load_block("blocks/A_result.bin",    comm, &A_A);
    load_block("blocks/A_B.bin",         comm, &A_B);

    load_block("blocks/B_neu_conj.bin",  comm, &B_neu);
    load_block("blocks/B_X_conj.bin",    comm, &B_X);
    load_block("blocks/B_A_conj.bin",    comm, &B_A);
    load_block("blocks/B_result.bin",    comm, &B_B);

    if (rank == 0)
        PetscPrintf(comm, "All blocks loaded.\n");

    Mat blocks[4][4] = {
        {neu_neu, neu_X,   neu_A,   neu_B},
        {X_neu,   X_X,     X_A,     X_B},
        {A_neu,   A_X,     A_A,     A_B},
        {B_neu,   B_X,     B_A,     B_B}
    };

    Mat Mnest;
    ierr = MatCreateNest(comm, 4, NULL, 4, NULL, (Mat*)blocks, &Mnest); CHKERRQ(ierr);
    ierr = MatNestSetVecType(Mnest, VECMPI); CHKERRQ(ierr);

    ierr = MatAssemblyBegin(Mnest, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Mnest, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    if (rank == 0)
        PetscPrintf(comm, "Converting nested matrix AIJ...\n");

    Mat A;
    ierr = MatConvert(Mnest, MATAIJ, MAT_INITIAL_MATRIX, &A); CHKERRQ(ierr);
    ierr = MatSetType(A, MATAIJ); CHKERRQ(ierr);

    if (rank == 0)
        PetscPrintf(comm, "Writing final AIJ matrix to: blocks/combined_matrix.bin\n");

    PetscViewer viewer;
    ierr = PetscViewerBinaryOpen(comm, "blocks/hamil.bin",
                                 FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    /*ierr = PetscViewerBinaryOpen(comm, "blocks/combined_matrix.bin",
                                 FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);*/
    ierr = MatView(A, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    if (rank == 0)
        PetscPrintf(comm, "DONE.\n");

    MatDestroy(&A);
    MatDestroy(&Mnest);
    MatDestroy(&neu_neu); MatDestroy(&neu_X);  MatDestroy(&neu_A);  MatDestroy(&neu_B);
    MatDestroy(&X_neu);   MatDestroy(&X_X);    MatDestroy(&X_A);    MatDestroy(&X_B);
    MatDestroy(&A_neu);   MatDestroy(&A_X);    MatDestroy(&A_A);    MatDestroy(&A_B);
    MatDestroy(&B_neu);   MatDestroy(&B_X);    MatDestroy(&B_A);    MatDestroy(&B_B);

    return 0;
}

