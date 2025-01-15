#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <omp.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <direct.h>

struct State {
    int s, e, n, l, j, m, t;
};

// Define a struct to store matrix elements
struct MatrixElement {
    int s1, s2, ss1, ss2;
    double me;
};

int main(int argc, char *argv[]) {
    // Check if the correct number of arguments is provided
    if (argc < 2) {
        printf("Usage: %s <Nmax>\n", argv[0]);
        return 1;
    }

    // Open the input file for reading
    FILE *file = fopen("../matrix_elements/hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_hwHO020.me2m", "r");
    if (!file) {
        printf("error: could not open file\n");
        return 1;
    }

    // Vector to store matrix elements
    std::vector<State> states;
    std::vector<MatrixElement> elements;

    // Variables to read element data
    int s, e, n, l, j, m, t;
    int s1, s2, ss1, ss2;
    double me;
    char line[256];

    // Skip all empty lines and comments until we reach the actual data
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '#' || line[0] == '\n') {
            continue;
        }
        if (sscanf(line, " %d %d %d %d %d %d %d", &s, &e, &n, &l, &j, &m, &t) == 7) {
            //printf("State: %d %d %d %d %d %d %d\n", s, e, n, l, j, m, t);
            states.push_back({s, e, n, l, j, m, t});
        } else if (sscanf(line, " %d %d %d %d %lf", &s1, &s2, &ss1, &ss2, &me) == 5) {
            //printf("Element: %d %d %d %d %lf\n", s1, s2, ss1, ss2, me);
            elements.push_back({s1, s2, ss1, ss2, me});
        }
    }

    // Close the input file
    fclose(file);

    // Check if states and elements are read correctly
    //printf("Number of states: %zu\n", states.size());
    //printf("Number of elements: %zu\n", elements.size());

    // Get the Nmax value from the command line arguments
    int Nmax = atoi(argv[1]);
    // Vector to store truncated elements
    std::vector<MatrixElement> truncated_elements;

    // Parallel loop to truncate elements based on Nmax
    int elements_size = elements.size();
    #pragma omp parallel for
    for (int i = 0; i < elements_size; i++) {
        if (states[elements[i].s1].e + states[elements[i].s2].e <= Nmax && states[elements[i].ss1].e + states[elements[i].ss2].e <= Nmax) {
            #pragma omp critical
            truncated_elements.push_back(elements[i]);
        }
    }

    // Debug print to check if elements are being truncated
    //printf("Number of truncated elements: %zu\n", truncated_elements.size());

    // for (const auto& element : truncated_elements) {
    //     printf("%d %d %d %d %lf\n", element.s1, element.s2, element.ss1, element.ss2, element.me);
    // }

    // Create the output directory if it doesn't exist
    const char* output_dir = "../matrix_elements/truncation";
    _mkdir(output_dir);

    // Construct the output file path
    char output_file_path[256];
    snprintf(output_file_path, sizeof(output_file_path), "%s/truncated_elements_Nmax%d.txt", output_dir, Nmax);
    
    // Open the output file for writing
    FILE *output_file = fopen(output_file_path, "w");
    if (!output_file) {
        printf("error: could not open output file\n");
        return 1;
    }

    // Write the truncated elements to the output file
    for (const auto& element : truncated_elements) {
        fprintf(output_file, "%d %d %d %d %lf\n", element.s1, element.s2, element.ss1, element.ss2, element.me);
    }

    // Close the output file
    fclose(output_file);

    return 0;
}
