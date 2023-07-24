#include "./include/includes.h"

/*-------------
Main function
--------------*/
int main(int argc, char *argv[])
{
    // Read parameters from command line
    if (argc != 5)
    {
        printf("Usage: %s <num_threads> <delta_t ODE (ms)> <delta_t PDE (ms)> <method>\n", argv[0]);
        exit(1);
    }

    // Get values from command line
    int numberThreads = atoi(argv[1]);
    float deltatODE = atof(argv[2]);
    float deltatPDE = atof(argv[3]);
    char *method = argv[4];

    // Call method
    char *options[] = {method, "true", "false", "false", "false", "false", "false"};
    runMethod(options, deltatODE, deltatPDE, numberThreads);

    return 0;
}
