#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "self_factors.h"

int read_states(const char *filename, State states[]) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening file");
        return 3;
    }

    int count = 0;
    char line[256];

    while (fgets(line, sizeof(line), fp) && count < MAX_STATES) {
        State s;
        char label[STATE_LABEL_LEN];
        int number;
        float energy;

        int matched = sscanf(line,
            "%9s %*[(]%*d%*[, ]%*d%*[, ]%*d%*[)] %d %f",
            label, &number, &energy);

        if (matched == 3) {
            strncpy(s.label, label, STATE_LABEL_LEN - 1);
            s.label[STATE_LABEL_LEN - 1] = '\0';
            s.number = number;
            s.energy = energy;
            states[count++] = s;
        }
    }

    fclose(fp);
    return count;
}

