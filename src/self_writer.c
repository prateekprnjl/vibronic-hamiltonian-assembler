/*Sparse writer*/

#include <petscmat.h>
#include <string.h>
#include "reader.h"
#include "self_factors.h"

PetscErrorCode build_self_matrix(const char *state_label,
                                 const State states[],
                                 int num_states,
                                 MPI_Comm comm)
{
    PetscErrorCode ierr;
    Mat baseMat, fullMat;

    int count = 0;
    for (int i = 0; i < num_states; i++) {
        if (strcmp(states[i].label, state_label) == 0) count++;
    }

    if (count == 0) {
        PetscPrintf(comm, "No states found for label %s\n", state_label);
        return 0;
    }

    const char *filename = NULL;
    if (strcmp(state_label, "neu") == 0) filename = "../data/blocks/Ham_blocks_block_0_0.bin";
    else if (strcmp(state_label, "X") == 0) filename = "../data/blocks/Ham_blocks_block_1_1.bin";
    else if (strcmp(state_label, "A") == 0) filename = "../data/blocks/Ham_blocks_block_2_2.bin";
    else if (strcmp(state_label, "B") == 0) filename = "../data/blocks/Ham_blocks_block_3_3.bin";
    else {
        PetscPrintf(comm, "Unknown state label: %s\n", state_label);
        return PETSC_ERR_USER;
    }

    baseMat = read_petsc_matrix(filename, comm);

    PetscInt nrows, ncols;
    ierr = MatGetSize(baseMat, &nrows, &ncols); CHKERRQ(ierr);

    PetscInt total_size = count * nrows;
    ierr = MatCreate(comm, &fullMat); CHKERRQ(ierr);
    ierr = MatSetSizes(fullMat, PETSC_DECIDE, PETSC_DECIDE, total_size, total_size); CHKERRQ(ierr);
    ierr = MatSetFromOptions(fullMat); CHKERRQ(ierr);
    ierr = MatSetUp(fullMat); CHKERRQ(ierr);

    PetscInt Istart, Iend;
    ierr = MatGetOwnershipRange(baseMat, &Istart, &Iend); CHKERRQ(ierr);

    int block_index = 0;
    for (int i = 0; i < num_states; i++) {
        if (strcmp(states[i].label, state_label) != 0) continue;

        PetscScalar energy_shift = states[i].energy * EV_TO_HARTREE;

        PetscInt row_offset = block_index * nrows;
        PetscInt col_offset = block_index * ncols;

        for (PetscInt r = Istart; r < Iend; r++) {
            const PetscInt *cols;
            const PetscScalar *vals;
            PetscInt ncols_row;

            ierr = MatGetRow(baseMat, r, &ncols_row, &cols, &vals); CHKERRQ(ierr);

            PetscInt *global_cols = malloc(ncols_row * sizeof(PetscInt));
            PetscScalar *new_vals = malloc(ncols_row * sizeof(PetscScalar));

            for (PetscInt c = 0; c < ncols_row; c++) {
                global_cols[c] = col_offset + cols[c];
                new_vals[c] = vals[c];
                if (cols[c] == r) new_vals[c] += energy_shift; 
            }

            PetscInt global_row = row_offset + r;
            ierr = MatSetValues(fullMat, 1, &global_row, ncols_row, global_cols, new_vals, INSERT_VALUES); CHKERRQ(ierr);

            free(global_cols);
            free(new_vals);
            ierr = MatRestoreRow(baseMat, r, &ncols_row, &cols, &vals); CHKERRQ(ierr);
        }

        block_index++;
    }

    ierr = MatAssemblyBegin(fullMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(fullMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    char outname[128];
    snprintf(outname, sizeof(outname), "blocks/%s_result.bin", state_label);
    PetscViewer viewer;
    ierr = PetscViewerBinaryOpen(comm, outname, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = MatView(fullMat, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    ierr = MatDestroy(&baseMat); CHKERRQ(ierr);
    ierr = MatDestroy(&fullMat); CHKERRQ(ierr);

    PetscPrintf(comm, "Wrote %s block matrix with %d blocks (%d-%d each %s)\n",
                state_label, count, (int)nrows, (int)ncols, outname);

    return 0;
}

