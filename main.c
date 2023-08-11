#include "./include/includes.h"

void resetSimulationParameters()
{
    // For stimulation
    Istim = 0.0;

    // For time step
    timeStepCounter = 0;
    timeStep = 0.0;

    // For execution time
    startTotal = 0.0;
    finishTotal = 0.0;
    elapsedTotal = 0.0;
    startODE = 0.0;
    finishODE = 0.0;
    elapsedODE = 0.0;
    startPDE = 0.0;
    finishPDE = 0.0;
    elapsedPDE = 0.0;

    // For velocity
    S1VelocityTag = true;
    S1Velocity = 0.0;
}

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
    resetSimulationParameters();
    runMethodMeasuringOnlyTotalTime(options, method, deltatODE, deltatPDE, numberThreads);

    // Measure Vulnerable Window
    // Update VWTag
    if (VWTag == false)
    {
        VWTag = true;
        stim2Begin = measureVWFrom;
        resetSimulationParameters();
        runMethod(options, method, deltatODE, deltatPDE, numberThreads);
    }
    // Update S2 begin
    while (stim2Begin + 1.0 <= measureVWTo)
    {
        stim2Begin += 1.0;
        resetSimulationParameters();
        runMethod(options, method, deltatODE, deltatPDE, numberThreads);
    }

    return 0;
}
