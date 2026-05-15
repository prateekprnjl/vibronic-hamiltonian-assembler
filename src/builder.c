#include <petscmat.h>
#include <string.h>
#include "fc_factors.h"

Mat read_petsc_matrix(const char *filename, MPI_Comm comm);

PetscErrorCode build_coupling_matrix(const char *state_label1,
                                     int base_state1,
                                     const char *state_label2,
                                     int base_state2,
                                     const char *base_filename,
                                     const FCEntry fc_entries[],
                                     int num_fc,
                                     MPI_Comm comm)
{
    PetscErrorCode ierr;
    Mat baseMat, fullMat;

    PetscPrintf(comm, "\nBuilding coupling matrix %s-%s from %s\n",
                state_label1, state_label2, base_filename);

    baseMat = read_petsc_matrix(base_filename, comm);
    PetscInt nrows, ncols;
    ierr = MatGetSize(baseMat, &nrows, &ncols); CHKERRQ(ierr);

    int max_num1 = -1, max_num2 = -1;
    for (int i = 0; i < num_fc; i++) {
        if (strcmp(fc_entries[i].state1, state_label1) == 0 &&
            strcmp(fc_entries[i].state2, state_label2) == 0)
        {
            if (fc_entries[i].num1 > max_num1) max_num1 = fc_entries[i].num1;
            if (fc_entries[i].num2 > max_num2) max_num2 = fc_entries[i].num2;
        }
    }

    if (max_num1 < 0 || max_num2 < 0) {
        PetscPrintf(comm, "No FC entries for %s-%s. Skipping.\n",
                    state_label1, state_label2);
        ierr = MatDestroy(&baseMat); CHKERRQ(ierr);
        return 0;
    }

    int count1 = max_num1 - base_state1 + 1;
    int count2 = max_num2 - base_state2 + 1;

    PetscInt total_rows = count1 * nrows;
    PetscInt total_cols = count2 * ncols;

    ierr = MatCreate(comm, &fullMat); CHKERRQ(ierr);
    ierr = MatSetSizes(fullMat, PETSC_DECIDE, PETSC_DECIDE, total_rows, total_cols); CHKERRQ(ierr);
    ierr = MatSetFromOptions(fullMat); CHKERRQ(ierr);
    ierr = MatSetUp(fullMat); CHKERRQ(ierr);

    PetscInt Istart, Iend;
    ierr = MatGetOwnershipRange(baseMat, &Istart, &Iend); CHKERRQ(ierr);

    for (int k = 0; k < num_fc; k++) {
        if (strcmp(fc_entries[k].state1, state_label1) != 0 ||
            strcmp(fc_entries[k].state2, state_label2) != 0)
            continue;

        int vib1 = fc_entries[k].num1 - base_state1;
        int vib2 = fc_entries[k].num2 - base_state2;
        double fc = fc_entries[k].fc_factor;

        PetscInt row_offset = vib1 * nrows;
        PetscInt col_offset = vib2 * ncols;

        for (PetscInt r = Istart; r < Iend; r++) {
            const PetscInt *cols;
            const PetscScalar *vals;
            PetscInt ncols_row;

            ierr = MatGetRow(baseMat, r, &ncols_row, &cols, &vals); CHKERRQ(ierr);

            PetscInt *global_cols = malloc(ncols_row * sizeof(PetscInt));
            PetscScalar *scaled_vals = malloc(ncols_row * sizeof(PetscScalar));

            for (PetscInt c = 0; c < ncols_row; c++) {
                global_cols[c] = col_offset + cols[c];
                scaled_vals[c] = vals[c] * fc;
            }

            PetscInt global_row = row_offset + r;
            ierr = MatSetValues(fullMat, 1, &global_row, ncols_row, global_cols, scaled_vals, ADD_VALUES); CHKERRQ(ierr);

            free(global_cols);
            free(scaled_vals);
            ierr = MatRestoreRow(baseMat, r, &ncols_row, &cols, &vals); CHKERRQ(ierr);
        }
    }

    ierr = MatAssemblyBegin(fullMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(fullMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    char outname[128];
    snprintf(outname, sizeof(outname), "blocks/%s_%s.bin", state_label1, state_label2);
    PetscViewer viewer;
    ierr = PetscViewerBinaryOpen(comm, outname, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = MatView(fullMat, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    Mat conjMat;
    ierr = MatHermitianTranspose(fullMat, MAT_INITIAL_MATRIX, &conjMat); CHKERRQ(ierr);

    snprintf(outname, sizeof(outname), "blocks/%s_%s_conj.bin", state_label2, state_label1);
    ierr = PetscViewerBinaryOpen(comm, outname, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = MatView(conjMat, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    ierr = MatDestroy(&conjMat); CHKERRQ(ierr);
    ierr = MatDestroy(&baseMat); CHKERRQ(ierr);
    ierr = MatDestroy(&fullMat); CHKERRQ(ierr);

    PetscPrintf(comm, "Wrote %s_%s and %s_%s_conj coupling matrices.\n",
                state_label1, state_label2, state_label2, state_label1);

    return 0;
}

