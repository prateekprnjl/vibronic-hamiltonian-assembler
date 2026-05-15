#include <stdio.h>
#include <petscmat.h>
#include "fc_factors.h"
#include "self_factors.h"
#include "reader.h"

int read_fc_data(const char *filename, FCEntry *data, int max_entries);
int read_states(const char *filename, State states[]);
Mat read_petsc_matrix(const char *filename, MPI_Comm comm);
PetscErrorCode assemble_full_matrix(MPI_Comm comm);
PetscErrorCode build_self_matrix(const char *state_label,
                                 const State states[],
                                 int num_states,
                                 MPI_Comm comm);
PetscErrorCode build_coupling_matrix(const char *state_label1,
                                     int base_state1,
                                     const char *state_label2,
                                     int base_state2,
                                     const char *base_filename,
                                     const FCEntry fc_entries[],
                                     int num_fc,
                                     MPI_Comm comm);

int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, NULL, NULL);

    State states[MAX_STATES];
    int num_states = read_states(INPUT_FILE2, states);

    build_self_matrix("X", states, num_states, PETSC_COMM_WORLD);
    build_self_matrix("A", states, num_states, PETSC_COMM_WORLD);
    build_self_matrix("B", states, num_states, PETSC_COMM_WORLD);
    build_self_matrix("neu", states, num_states, PETSC_COMM_WORLD);

    FCEntry overlaps[MAX_ENTRIES];
    int n_overlaps = read_fc_data(INPUT_FILE, overlaps, MAX_ENTRIES);

    MPI_Comm comm = PETSC_COMM_WORLD;

    build_coupling_matrix("neu", 0, "X", 1, "../data/blocks/Ham_blocks_block_0_1.bin", overlaps, n_overlaps, comm);
    build_coupling_matrix("neu", 0, "A", 15, "../data/blocks/Ham_blocks_block_0_2.bin", overlaps, n_overlaps, comm);
    build_coupling_matrix("neu", 0, "B", 29, "../data/blocks/Ham_blocks_block_0_3.bin", overlaps, n_overlaps, comm);
    build_coupling_matrix("X", 1, "A", 15, "../data/blocks/Ham_blocks_block_1_2.bin", overlaps, n_overlaps, comm);
    build_coupling_matrix("X", 1, "B", 29, "../data/blocks/Ham_blocks_block_1_3.bin", overlaps, n_overlaps, comm);
    build_coupling_matrix("A", 15, "B", 29, "../data/blocks/Ham_blocks_block_2_3.bin", overlaps, n_overlaps, comm);

    assemble_full_matrix(comm);

    PetscFinalize();

    return 0;
}

