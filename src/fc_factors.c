#include "fc_factors.h"

int read_fc_data(const char *filename, FCEntry *data, int max_entries) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening file");
        exit(2);
    }

    int count = 0;
    char line[512];

    while (fgets(line, sizeof(line), fp) && count < max_entries) {
        FCEntry e;
        memset(&e, 0, sizeof(FCEntry));

        int matched = sscanf(line,
            "%7s %*[(]%*d%*[, ]%*d%*[, ]%*d%*[)] %7s %*[(]%*d%*[, ]%*d%*[, ]%*d%*[)] %d %d [%lf %lf] %*lf",
            e.state1, e.state2, &e.num1, &e.num2,
            &e.energy1, &e.fc_factor);

        if (matched == 6)
            data[count++] = e;
    }
    fclose(fp);
    return (count);
}

