# Vibronic Hamiltonian assembler

A C-based scientific computing program for assembling diabatic/vibronic Hamiltonian matrices (or Dipoles with small modifications) from precomputed Hamiltonian blocks using Franck–Condon (FC) factor approximations.

---

## Overview

This program constructs large Hamiltonian matrices by combining smaller block Hamiltonians with coupling terms derived from Franck–Condon factors.

The workflow is intended for vibronic coupling and multi-state Hamiltonian assembly calculations commonly used in quantum dynamics, molecular spectroscopy, and nonadiabatic simulations.

---

## Features

- Reads binary Hamiltonian block matrices
- Applies Franck–Condon factor approximations
- Constructs coupled Hamiltonian matrices
- Supports self-coupling and off-diagonal coupling terms
- Modular C implementation for easy extension

---

## Project Structure

```text
src/
    assemble.c         Matrix assembly routines
    builder.c          Hamiltonian construction logic
    fc_factors.c       Franck–Condon factor handling
    self_factors.c     Self-coupling factor generation
    self_writer.c      Output writer utilities
    reader.c           Input readers
    main.c             Main driver

data/
    blocks/            Hamiltonian block files
    fc/                FC factor input files

docs/
    theory.pdf          Theory for the Hamiltonian assembler. (logic for dipoles as well, present program can calculate that with minor modifications)
```

---

## Input Files

### Hamiltonian Blocks

Example:

```text
Ham_blocks_block_0_0.bin
Ham_blocks_block_0_1.bin
...
```

These binary files contain precomputed Hamiltonian matrix blocks.

This repository does not include the full production Hamiltonian block datasets due to their large size. Users are expected to generate or supply their own Hamiltonian block data and Franck–Condon factor files compatible with the formats described in this repository.

### FC Factor Files

Example:

```text
fc_factor.out
self_fc_factor.out
```

These files contain Franck–Condon overlap/coupling factors used during assembly.

---

## Compilation

```bash
make
```

or manually:

```bash
gcc -O2 src/*.c -o ham_builder
```

---

## Usage

Example:

```bash
./ham_builder
```

Modify input paths and parameters directly in `main.c` or configuration sections as needed.

---

## Theory

The Hamiltonian is assembled using blockwise coupling expressions weighted by Franck–Condon overlap factors.

A schematic form is:

$H_{ij}^{total} = H_{ij}^{block} \times FC_{ij}$

where:

- $H_{ij}^{block}$ is the electronic/vibronic block contribution
- $FC_{ij}$ is the Franck–Condon overlap factor

Additional self-coupling and diagonal corrections may also be included. Further details are provided in docs/theory.md

---

## Output

The program generates assembled Hamiltonian matrices suitable for downstream quantum dynamics or spectroscopy calculations.

Output formats may include:

- Binary matrices
- Intermediate coupling files
- Diagnostic logs

---

## Dependencies

- GCC or compatible C compiler
- Standard C libraries
- PETSc, SLEPc

---

## Future Improvements

- Configuration file support
- Python interface
- HDF5-native workflows

---
## License

MIT License
