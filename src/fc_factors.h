#ifndef FC_FACTORS_H
#define FC_FACTORS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ENTRIES 2000
#define MAX_LABEL_LEN 8
#define INPUT_FILE "../data/fc/fc_factor.out"
#define EV_TO_HARTREE 0.0367493

typedef struct {
    char state1[MAX_LABEL_LEN];
    char state2[MAX_LABEL_LEN];
    int num1;
    int num2;
    double energy1;
    double fc_factor;
} FCEntry;

int read_fc_data(const char *filename, FCEntry *data, int max_entries);

#endif

