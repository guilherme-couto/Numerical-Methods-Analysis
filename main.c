#include "./include/includes.h"

/*-------------
Main function
--------------*/
int main(int argc, char *argv[])
{
    // Read parameters from command line
    if (argc != 6)
    {
        printf("Usage: %s <num_threads> <delta_t ODE (ms)> <delta_t PDE (ms)> <method> <fibrosis (0-false or 1-true)>\n", argv[0]);
        exit(1);
    }

    // Get values from command line
    int numberThreads = atoi(argv[1]);
    float deltatODE = atof(argv[2]);
    float deltatPDE = atof(argv[3]);
    char *method = argv[4];
    bool fibrosis = atoi(argv[5]);

    // Call method
    if (!fibrosis)
    {
        fibrosisFactor = 1.0;
    }
    bool options[] = {fibrosis, true, false, false, false, false, false};
    runMethod(options, method, deltatODE, deltatPDE, numberThreads);

    return 0;
}
