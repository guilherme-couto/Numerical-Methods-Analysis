#include "./include/includes.h"

/*-------------
Main function
--------------*/
int main(int argc, char *argv[])
{
    // Read parameters from command line
    if (argc != 8)
    {
        // For convergence analysis
        if (argc == 7)
        {
            int numberThreads = atoi(argv[1]);
            float delta_x = atof(argv[2]);
            float delta_t = atof(argv[3]);
            char *method = argv[4];
            int rate = atoi(argv[5]);
            int spatialRate = atoi(argv[6]);

            // Call method
            bool fibrosis = false;
            if (!fibrosis)
            {
                fibrosisFactor = 1.0;
            }
            //runMethodConvergenceAnalysis(fibrosis, method, delta_x, delta_t, numberThreads, spatialRate, rate);
            return 0;
        }
        printf("Usage: %s <num_threads> <delta_t ODE (ms)> <delta_t PDE (ms)> <method> <fibrosis (0-false or 1-true)> <vulnerability window (0-false or 1-true)> <only speed (0-false or 1-true)>\n", argv[0]);
        exit(1);
    }

    // Get values from command line
    int numberThreads = atoi(argv[1]);
    float deltatODE = atof(argv[2]);
    float deltatPDE = atof(argv[3]);
    char *method = argv[4];
    bool fibrosis = atoi(argv[5]);
    bool VWmeasure = atoi(argv[6]);
    bool onlySpeed = atoi(argv[7]);

    // Call method
    if (!fibrosis)
    {
        fibrosisFactor = 1.0;
    }

    // bool haveFibrosis = options[0];
    // bool measureTotalTime = options[1];
    // bool saveDataToError = options[2];
    // bool saveDataToGif = options[3];
    // bool measureTimeParts = options[4];
    // bool measureS1Velocity = options[5];
    bool options[] = {fibrosis, true, true, true, true, true};

    if (onlySpeed == false)
    {
        runMethod(options, method, deltatODE, deltatPDE, numberThreads);
    }
    else if (onlySpeed == true)
    {
        options[1] = false;
        options[2] = false;
        options[3] = false;
        runMethod(options, method, deltatODE, deltatPDE, numberThreads);
        resetSimulationParameters();
        runMethodOnlyTotalTimeAnalysis(options, method, deltatODE, deltatPDE, numberThreads);
    }
    

    // Measure Vulnerable Window
    if (VWmeasure == true)
    {
        // bool haveFibrosis = options[0];
        // bool measureTotalTime = options[1];
        // bool saveDataToError = options[2];
        // bool saveDataToGif = options[3];
        // bool measureTimeParts = options[4];
        // bool measureS1Velocity = options[5];
        bool options[] = {fibrosis, true, false, false, true, true};

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
    }

    return 0;
}
