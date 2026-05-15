#ifndef SELF_READER_H
#define SELF_READER_H

#define MAX_STATES 100
#define STATE_LABEL_LEN 10
#define EV_TO_HARTREE 0.0367493

#define INPUT_FILE2 "../data/fc/self_fc_factor.out"

typedef struct {
    char label[STATE_LABEL_LEN];
    int number;
    float energy;
} State;

int read_states(const char *filename, State states[]);

#endif

