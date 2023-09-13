#ifndef METHODS_H
#define METHODS_H

#include "includes.h"

// For stimulation
double Istim = 0.0;

// For time step
int timeStepCounter = 0;
double timeStep = 0.0;

// For execution time
double startTotal = 0.0, finishTotal = 0.0, elapsedTotal = 0.0;
double startODE = 0.0, finishODE = 0.0, elapsedODE = 0.0;
double startPDE = 0.0, finishPDE = 0.0, elapsedPDE = 0.0;
double lastCheckpointTime = 0.0;

// For velocity
bool S1VelocityTag = true;
double S1Velocity = 0.0;

// For vulnerability window
bool VWTag = false;
double measureVWFrom = 117.0, measureVWTo = 132.0;
double lowerVWBound = 999.0, upperVWBound = 0.0;


//############################################
//##                                        ##
//##     Adapted FitzHugh-Nagumo (AFHN)     ##
//##                                        ##
//############################################
#if defined(AFHN)
void runMethod(bool options[], char *method, float deltatODE, float deltatPDE, int numberThreads)
{
    // Get options
    bool haveFibrosis = options[0];
    bool measureTotalTime = options[1];
    bool saveDataToError = options[2];
    bool saveDataToGif = options[3];
    bool measureTimeParts = options[4];
    bool measureS1Velocity = options[5];

    // Number of steps
    int N = round(L / deltax) + 1;                          // Spatial steps (square tissue)
    int M = round(T / deltatPDE) + 1;                       // Number of time steps
    int PdeOdeRatio = round(deltatPDE / deltatODE);     // Ratio between PDE and ODE time steps

    // Allocate and populate time array
    double *time;
    // allocateAndPopulateTime(M, time, deltatPDE);
    time = (double *)malloc(M * sizeof(double));
    for (int i = 0; i < M; i++)
    {
        time[i] = i * deltatPDE;
    }

    // Allocate and initialize variables
    double **V, **W;
    V = (double **)malloc(N * sizeof(double *));
    W = (double **)malloc(N * sizeof(double *));

    for (int i = 0; i < N; i++)
    {
        V[i] = (double *)malloc(N * sizeof(double));
        W[i] = (double *)malloc(N * sizeof(double));
    }
    initializeVariables(N, V, W);

    // Diffusion coefficient - isotropic
    double D = sigma / (chi * Cm);
    double phi = D * deltatPDE / (deltax * deltax);    // For Thomas algorithm - isotropic

    // Variables
    int i, j;                                          // i for y-axis and j for x-axis
    double actualV, actualW;
    double Rw;
    double **Vtilde, **Wtilde, **Rv, **rightside, **solution;
    Vtilde = (double **)malloc(N * sizeof(double *));
    Wtilde = (double **)malloc(N * sizeof(double *));
    Rv = (double **)malloc(N * sizeof(double *));
    rightside = (double **)malloc(N * sizeof(double *));
    solution = (double **)malloc(N * sizeof(double *));

    // Auxiliary arrays for Thomas algorithm 2nd order approximation
    double **c_ = (double **)malloc((N) * sizeof(double *));
    double **d_ = (double **)malloc((N) * sizeof(double *));
    for (int i = 0; i < N; i++)
    {
        Vtilde[i] = (double *)malloc(N * sizeof(double));
        Wtilde[i] = (double *)malloc(N * sizeof(double));
        Rv[i] = (double *)malloc(N * sizeof(double));
        rightside[i] = (double *)malloc(N * sizeof(double));
        solution[i] = (double *)malloc(N * sizeof(double));
        c_[i] = (double *)malloc(N * sizeof(double));
        d_[i] = (double *)malloc(N * sizeof(double));
    }

    // Discretized limits of stimulation area
    int discS1xLimit = round(stim1xLimit / deltax);
    int discS1yLimit = round(stim1yLimit / deltay);
    int discS2xMax = round(stim2xMax / deltax);
    int discS2xMin = round(stim2xMin / deltax);
    int discS2yMax = N;
    int discS2yMin = N - round(stim2yMax / deltay);

    // Discritized limits of fibrotic area
    int discFibxMax = round(fibrosisMaxX / deltax);
    int discFibxMin = round(fibrosisMinX / deltax);
    int discFibyMax = N - round(fibrosisMinY / deltay);
    int discFibyMin = N - round(fibrosisMaxY / deltay);

    // File names
    char framesFileName[MAX_STRING_SIZE], infosFileName[MAX_STRING_SIZE];
    sprintf(framesFileName, "frames-%d-%.3lf-%.3lf.txt", numberThreads, deltatODE, deltatPDE);
    sprintf(infosFileName, "infos-%d-%.3lf-%.3lf.txt", numberThreads, deltatODE, deltatPDE);
    int saverate = ceil(M / 100.0);
    FILE *fpFrames, *fpInfos;

    // Create directories and files
    char pathToSaveData[MAX_STRING_SIZE];
    if (haveFibrosis)
    {
        createDirectories(pathToSaveData, method, "AFHN-Fibro");
    }
    else
    {
        createDirectories(pathToSaveData, method, "AFHN");
    }
    
    // File pointers
    char aux[MAX_STRING_SIZE];
    if (VWTag == false)
    {
        
        sprintf(aux, "%s/%s", pathToSaveData, infosFileName);
        fpInfos = fopen(aux, "w");
    }
    else
    {
        sprintf(aux, "%s/%s", pathToSaveData, infosFileName);
        fpInfos = fopen(aux, "a");
    }
    if (saveDataToGif == false)
    {
        sprintf(aux, "%s/%s", pathToSaveData, framesFileName);
        fpFrames = fopen(aux, "a");
    }
    else
    {
        sprintf(aux, "%s/%s", pathToSaveData, framesFileName);
        fpFrames = fopen(aux, "w");
    }
    
    /*--------------------
    --  ADI 1st order   --
    ----------------------*/
    if (strcmp(method, "ADI1") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, S1Velocity, S1VelocityTag, VWTag, saverate, fpFrames, saveDataToGif, \
        Rv, rightside, solution, startODE, finishODE, elapsedODE, startPDE, finishPDE, elapsedPDE)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Start measuring ODE execution time
                #pragma omp master
                {
                    startODE = omp_get_wtime();
                }

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V and W without diffusion
                            V[i][j] = actualV + deltatODE * (reactionV(actualV, actualW) + Istim);
                            W[i][j] = actualW + deltatODE * reactionW(actualV, actualW);

                            // Update right side of Thomas algorithm
                            rightside[j][i] = V[i][j];
                        }
                    }
                }

                // Finish measuring ODE execution time and start measuring PDE execution time
                #pragma omp master
                {
                    finishODE = omp_get_wtime();
                    elapsedODE += finishODE - startODE;
                    startPDE = omp_get_wtime();
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines)
                #pragma omp barrier
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, phi, c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, phi, c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        Vtilde[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns)
                #pragma omp barrier
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(Vtilde[i], V[i], N, phi, c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(Vtilde[i], V[i], N, phi, c_[i], d_[i], N, 0);
                    }
                }

                // Finish measuring PDE execution time
                #pragma omp master
                {
                    finishPDE = omp_get_wtime();
                    elapsedPDE += finishPDE - startPDE;
                }

                // Save frames
                #pragma omp master
                {
                    if (VWTag == false)
                    {
                        // Write frames to file
                        if (timeStepCounter % saverate == 0 && saveDataToGif == true)
                        {
                            fprintf(fpFrames, "%lf\n", time[timeStepCounter]);
                            for (i = 0; i < N; i++)
                            {
                                for (j = 0; j < N; j++)
                                {
                                    fprintf(fpFrames, "%lf ", V[i][j]);
                                }
                                fprintf(fpFrames, "\n");
                            }
                        }

                        // Check S1 velocity
                        if (S1VelocityTag)
                        {
                            if (V[0][N-1] >= 80)
                            {
                                S1Velocity = ((10 * (L - stim1xLimit)) / (time[timeStepCounter]));
                                S1VelocityTag = false;
                            }
                        }
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }

    
    /*-------------------------
    --  SSI-ADI (2nd order)  --
    --------------------------*/
    else if (strcmp(method, "SSI-ADI") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW, Rw) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, S1Velocity, S1VelocityTag, VWTag, saverate, fpFrames, saveDataToGif, \
        Rv, rightside, solution, startODE, finishODE, elapsedODE, startPDE, finishPDE, elapsedPDE)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Start measuring ODE execution time
                #pragma omp master
                {
                    startODE = omp_get_wtime();
                }

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion and W without diffusion (Heun)
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));
                            
                            Rw = deltatODE * reactionW(actualV, actualW);
                            //Wtilde[i][j] = actualW + (0.5 * deltatODE * reactionW(actualV, actualW));
                            Wtilde[i][j] = actualW + 0.5 * Rw;

                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], Wtilde[i][j]) + Istim);

                            // Update W explicitly (Heun)
                            W[i][j] = actualW + (deltatODE * reactionW(Vtilde[i][j], Wtilde[i][j]));
                        }
                    }
                }

                // Finish measuring ODE execution time and start measuring PDE execution time
                #pragma omp master
                {
                    finishODE = omp_get_wtime();
                    elapsedODE += finishODE - startODE;
                    startPDE = omp_get_wtime();
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Finish measuring PDE execution time
                #pragma omp master
                {
                    finishPDE = omp_get_wtime();
                    elapsedPDE += finishPDE - startPDE;
                }

                // Save frames
                #pragma omp master
                {
                    if (VWTag == false)
                    {
                        // Write frames to file
                        if (timeStepCounter % saverate == 0 && saveDataToGif == true)
                        {
                            fprintf(fpFrames, "%lf\n", time[timeStepCounter]);
                            for (i = 0; i < N; i++)
                            {
                                for (j = 0; j < N; j++)
                                {
                                    fprintf(fpFrames, "%lf ", V[i][j]);
                                }
                                fprintf(fpFrames, "\n");
                            }
                        }

                        // Check S1 velocity
                        if (S1VelocityTag)
                        {
                            if (V[0][N-1] >= 80)
                            {
                                S1Velocity = ((10 * (L - stim1xLimit)) / (time[timeStepCounter]));
                                S1VelocityTag = false;
                            }
                        }
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    /*-------------------------
    --  ADI1.5 (1.5 order)  --
    --------------------------*/
    else if (strcmp(method, "ADI1.5") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, S1Velocity, S1VelocityTag, VWTag, saverate, fpFrames, saveDataToGif, \
        Rv, rightside, solution, startODE, finishODE, elapsedODE, startPDE, finishPDE, elapsedPDE)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Start measuring ODE execution time
                #pragma omp master
                {
                    startODE = omp_get_wtime();
                }

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));

                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], actualW) + Istim);

                            // Update W explicitly
                            W[i][j] = actualW + (deltatODE * reactionW(Vtilde[i][j], actualW));
                        }
                    }
                }

                // Finish measuring ODE execution time and start measuring PDE execution time
                #pragma omp master
                {
                    finishODE = omp_get_wtime();
                    elapsedODE += finishODE - startODE;
                    startPDE = omp_get_wtime();
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                        
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Finish measuring PDE execution time
                #pragma omp master
                {
                    finishPDE = omp_get_wtime();
                    elapsedPDE += finishPDE - startPDE;
                }

                // Save frames
                #pragma omp master
                {
                    if (VWTag == false)
                    {
                        // Write frames to file
                        if (timeStepCounter % saverate == 0 && saveDataToGif == true)
                        {
                            fprintf(fpFrames, "%lf\n", time[timeStepCounter]);
                            for (i = 0; i < N; i++)
                            {
                                for (j = 0; j < N; j++)
                                {
                                    fprintf(fpFrames, "%lf ", V[i][j]);
                                }
                                fprintf(fpFrames, "\n");
                            }
                        }

                        // Check S1 velocity
                        if (S1VelocityTag)
                        {
                            if (V[0][N-1] >= 80)
                            {
                                S1Velocity = ((10 * (L - stim1xLimit)) / (time[timeStepCounter]));
                                S1VelocityTag = false;
                            }
                        }
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    /*-------------------------
    --  ADI1.5 (1.5.2 order)  --
    --------------------------*/
    else if (strcmp(method, "ADI1.5.2") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, S1Velocity, S1VelocityTag, VWTag, saverate, fpFrames, saveDataToGif, \
        Rv, rightside, solution, startODE, finishODE, elapsedODE, startPDE, finishPDE, elapsedPDE)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Start measuring ODE execution time
                #pragma omp master
                {
                    startODE = omp_get_wtime();
                }

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));

                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], actualW) + Istim);

                            // Update W explicitly
                            W[i][j] = actualW + (deltatODE * reactionW(actualV, actualW));
                        }
                    }
                }

                // Finish measuring ODE execution time and start measuring PDE execution time
                #pragma omp master
                {
                    finishODE = omp_get_wtime();
                    elapsedODE += finishODE - startODE;
                    startPDE = omp_get_wtime();
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                        
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Finish measuring PDE execution time
                #pragma omp master
                {
                    finishPDE = omp_get_wtime();
                    elapsedPDE += finishPDE - startPDE;
                }

                // Save frames
                #pragma omp master
                {
                    if (VWTag == false)
                    {
                        // Write frames to file
                        if (timeStepCounter % saverate == 0 && saveDataToGif == true)
                        {
                            fprintf(fpFrames, "%lf\n", time[timeStepCounter]);
                            for (i = 0; i < N; i++)
                            {
                                for (j = 0; j < N; j++)
                                {
                                    fprintf(fpFrames, "%lf ", V[i][j]);
                                }
                                fprintf(fpFrames, "\n");
                            }
                        }

                        // Check S1 velocity
                        if (S1VelocityTag)
                        {
                            if (V[0][N-1] >= 80)
                            {
                                S1Velocity = ((10 * (L - stim1xLimit)) / (time[timeStepCounter]));
                                S1VelocityTag = false;
                            }
                        }
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    /*-------------------------
    --  ADI1.5 (1.5.3 order)  --
    --------------------------*/
    else if (strcmp(method, "ADI1.5.3") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW, Rw) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, S1Velocity, S1VelocityTag, VWTag, saverate, fpFrames, saveDataToGif, \
        Rv, rightside, solution, startODE, finishODE, elapsedODE, startPDE, finishPDE, elapsedPDE)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Start measuring ODE execution time
                #pragma omp master
                {
                    startODE = omp_get_wtime();
                }

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));
                            
                            Rw = deltatODE * reactionW(actualV, actualW);
                            Wtilde[i][j] = actualW + 0.5 * Rw;
                            
                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], Wtilde[i][j]) + Istim);

                            // Update W explicitly
                            W[i][j] = actualW + Rw;
                            // W[i][j] = actualW + (deltatODE * reactionW(Vtilde[i][j], Wtilde[i][j]));
                        }
                    }
                }

                // Finish measuring ODE execution time and start measuring PDE execution time
                #pragma omp master
                {
                    finishODE = omp_get_wtime();
                    elapsedODE += finishODE - startODE;
                    startPDE = omp_get_wtime();
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Finish measuring PDE execution time
                #pragma omp master
                {
                    finishPDE = omp_get_wtime();
                    elapsedPDE += finishPDE - startPDE;
                }

                // Save frames
                #pragma omp master
                {
                    if (VWTag == false)
                    {
                        // Write frames to file
                        if (timeStepCounter % saverate == 0 && saveDataToGif == true)
                        {
                            fprintf(fpFrames, "%lf\n", time[timeStepCounter]);
                            for (i = 0; i < N; i++)
                            {
                                for (j = 0; j < N; j++)
                                {
                                    fprintf(fpFrames, "%lf ", V[i][j]);
                                }
                                fprintf(fpFrames, "\n");
                            }
                        }

                        // Check S1 velocity
                        if (S1VelocityTag)
                        {
                            if (V[0][N-1] >= 80)
                            {
                                S1Velocity = ((10 * (L - stim1xLimit)) / (time[timeStepCounter]));
                                S1VelocityTag = false;
                            }
                        }
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }



    /*--------------------
    --  Forward Euler   --
    ----------------------*/
    else if (strcmp(method, "FE") == 0)
    {
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, S1Velocity, S1VelocityTag, VWTag, saverate, fpFrames, saveDataToGif, \
        Rv, rightside, solution, startODE, finishODE, elapsedODE, startPDE, finishPDE, elapsedPDE)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Start measuring ODE execution time
                #pragma omp master
                {
                    startODE = omp_get_wtime();
                }

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];

                            // Update V and W
                            Vtilde[i][j] = actualV + deltatODE * (reactionV(actualV, actualW) + Istim);
                            V[i][j] = Vtilde[i][j];
                            W[i][j] = actualW + deltatODE * reactionW(actualV, actualW);
                        }
                    }
                }

                // Finish measuring ODE execution time and start measuring PDE execution time
                #pragma omp master
                {
                    finishODE = omp_get_wtime();
                    elapsedODE += finishODE - startODE;
                    startPDE = omp_get_wtime();
                }

                // Resolve PDEs (Diffusion)
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        V[i][j] = Vtilde[i][j] + phi * (iDiffusion2nd(i, j, N, Vtilde, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, Vtilde, discFibxMax, discFibxMin, discFibyMax, discFibyMin));
                    }
                }

                // Finish measuring PDE execution time
                #pragma omp master
                {
                    finishPDE = omp_get_wtime();
                    elapsedPDE += finishPDE - startPDE;
                }

                // Save frames
                #pragma omp master
                {
                    if (VWTag == false)
                    {
                        // Write frames to file
                        if (timeStepCounter % saverate == 0 && saveDataToGif == true)
                        {
                            fprintf(fpFrames, "%lf\n", time[timeStepCounter]);
                            for (i = 0; i < N; i++)
                            {
                                for (j = 0; j < N; j++)
                                {
                                    fprintf(fpFrames, "%lf ", V[i][j]);
                                }
                                fprintf(fpFrames, "\n");
                            }
                        }

                        // Check S1 velocity
                        if (S1VelocityTag)
                        {
                            if (V[0][N-1] >= 80)
                            {
                                S1Velocity = ((10 * (L - stim1xLimit)) / (time[timeStepCounter]));
                                S1VelocityTag = false;
                            }
                        }
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    // Write infos to file
    if (VWTag == false)
    {   
        fprintf(fpInfos, "S1 velocity: %lf m/s\n", S1Velocity);
        fprintf(fpInfos, "S2 begin: %.1lf ms\n", stim2Begin);
        fprintf(fpInfos, "PDE/ODE ratio: %d\n", PdeOdeRatio);
        fprintf(fpInfos, "ODE execution time: %lf seconds\n", elapsedODE);
        fprintf(fpInfos, "PDE execution time: %lf seconds\n", elapsedPDE);
        if (haveFibrosis)
        {
            fprintf(fpInfos, "Fibrosis factor: %.2lf\n", fibrosisFactor);
            fprintf(fpInfos, "Fibrosis region: (%.2lf, %.2lf) to (%.2lf, %.2lf)\n", fibrosisMinX, fibrosisMinY, fibrosisMaxX, fibrosisMaxY);
        }

        if (saveDataToError == true)
        {
            char lastFrameFileName[MAX_STRING_SIZE];
            sprintf(lastFrameFileName, "last-%d-%.3lf-%.3lf.txt", numberThreads, deltatODE, deltatPDE);
            FILE *fpLast;
            sprintf(aux, "%s/%s", pathToSaveData, lastFrameFileName);
            fpLast = fopen(aux, "w");
            FILE *fpSpiralV = fopen("lastSpiralV.txt", "w");
            FILE *fpSpiralW = fopen("lastSpiralW.txt", "w");
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    fprintf(fpLast, "%lf ", V[i][j]);
                    fprintf(fpSpiralV, "%lf ", V[i][j]);
                    fprintf(fpSpiralW, "%lf ", W[i][j]);
                }
                fprintf(fpLast, "\n");
                fprintf(fpSpiralV, "\n");
                fprintf(fpSpiralW, "\n");
            }
            fclose(fpLast);   
            fclose(fpSpiralV);
            fclose(fpSpiralW);  
        }   
    }

    // If vulnerability window is being measure, write last frame to file
    if (VWTag == true)
    {
        bool auxTag = true;
        bool above80 = false;

        // First, check if there is any voltage in the last frame with negative or nan value
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++) 
            {
                if (V[i][j] < 0 || isnan(V[i][j]) != 0)
                {
                    auxTag = false;
                    break;
                }
                if (V[i][j] >= 80.0)
                {
                    above80 = true;
                }
            }
            if (auxTag == false)
            {
                break;
            }
        }

        if (auxTag && above80)
        {
            // Update VW boundaries
            if (stim2Begin < lowerVWBound)
            {
                lowerVWBound = stim2Begin;
            }
            if (stim2Begin > upperVWBound)
            {
                upperVWBound = stim2Begin;
            }
        }

        if (stim2Begin + 1.0 > measureVWTo)
        {
            // Write infos to file
            fprintf(fpInfos, "Vulnerability Window: [%.1lf, %.1lf] ms\n", lowerVWBound, upperVWBound);
        }
    }

    // Close files
    fclose(fpFrames);
    fclose(fpInfos);

    // Free memory
    free(time);
    free(V);
    free(W);
    free(Vtilde);
    free(Wtilde);
    free(Rv);
    free(rightside);
    free(solution);
    free(c_);
    free(d_);
}

void runMethodOnlyTotalTimeAnalysis(bool options[], char *method, float deltatODE, float deltatPDE, int numberThreads)
{
    // Get options
    bool haveFibrosis = options[0];

    // Number of steps
    int N = round(L / deltax) + 1;                          // Spatial steps (square tissue)
    int M = round(T / deltatPDE) + 1;                       // Number of time steps
    int PdeOdeRatio = round(deltatPDE / deltatODE);     // Ratio between PDE and ODE time steps

    // Allocate and populate time array
    double *time;
    // allocateAndPopulateTime(M, time, deltatPDE);
    time = (double *)malloc(M * sizeof(double));
    for (int i = 0; i < M; i++)
    {
        time[i] = i * deltatPDE;
    }

    // Allocate and initialize variables
    double **V, **W;
    V = (double **)malloc(N * sizeof(double *));
    W = (double **)malloc(N * sizeof(double *));

    for (int i = 0; i < N; i++)
    {
        V[i] = (double *)malloc(N * sizeof(double));
        W[i] = (double *)malloc(N * sizeof(double));
    }
    initializeVariables(N, V, W);

    // Diffusion coefficient - isotropic
    double D = sigma / (chi * Cm);
    double phi = D * deltatPDE / (deltax * deltax);    // For Thomas algorithm - isotropic

    // Variables
    int i, j;                                          // i for y-axis and j for x-axis
    double actualV, actualW;
    double Rw;
    double **Vtilde, **Wtilde, **Rv, **rightside, **solution;
    Vtilde = (double **)malloc(N * sizeof(double *));
    Wtilde = (double **)malloc(N * sizeof(double *));
    Rv = (double **)malloc(N * sizeof(double *));
    rightside = (double **)malloc(N * sizeof(double *));
    solution = (double **)malloc(N * sizeof(double *));

    // Auxiliary arrays for Thomas algorithm 2nd order approximation
    double **c_ = (double **)malloc((N) * sizeof(double *));
    double **d_ = (double **)malloc((N) * sizeof(double *));
    for (int i = 0; i < N; i++)
    {
        Vtilde[i] = (double *)malloc(N * sizeof(double));
        Wtilde[i] = (double *)malloc(N * sizeof(double));
        Rv[i] = (double *)malloc(N * sizeof(double));
        rightside[i] = (double *)malloc(N * sizeof(double));
        solution[i] = (double *)malloc(N * sizeof(double));
        c_[i] = (double *)malloc(N * sizeof(double));
        d_[i] = (double *)malloc(N * sizeof(double));
    }

    // Discretized limits of stimulation area
    int discS1xLimit = round(stim1xLimit / deltax);
    int discS1yLimit = round(stim1yLimit / deltay);
    int discS2xMax = round(stim2xMax / deltax);
    int discS2xMin = round(stim2xMin / deltax);
    int discS2yMax = N;
    int discS2yMin = N - round(stim2yMax / deltay);

    // Discritized limits of fibrotic area
    int discFibxMax = round(fibrosisMaxX / deltax);
    int discFibxMin = round(fibrosisMinX / deltax);
    int discFibyMax = N - round(fibrosisMinY / deltay);
    int discFibyMin = N - round(fibrosisMaxY / deltay);

    // Create directories and files
    char pathToSaveData[MAX_STRING_SIZE];
    if (haveFibrosis)
    {
        createDirectories(pathToSaveData, method, "AFHN-Fibro");
    }
    else
    {
        createDirectories(pathToSaveData, method, "AFHN");
    }


    /*--------------------
    --  ADI 1st order   --
    ----------------------*/
    if (strcmp(method, "ADI1") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V and W without diffusion
                            V[i][j] = actualV + deltatODE * (reactionV(actualV, actualW) + Istim);
                            W[i][j] = actualW + deltatODE * reactionW(actualV, actualW);

                            // Update right side of Thomas algorithm
                            rightside[j][i] = V[i][j];
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines)
                #pragma omp barrier
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, phi, c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, phi, c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        Vtilde[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns)
                #pragma omp barrier
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(Vtilde[i], V[i], N, phi, c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(Vtilde[i], V[i], N, phi, c_[i], d_[i], N, 0);
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }

    
    /*-------------------------
    --  SSI-ADI (2nd order)  --
    --------------------------*/
    else if (strcmp(method, "SSI-ADI") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion and W without diffusion (Heun)
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));
                            Wtilde[i][j] = actualW + (0.5 * deltatODE * reactionW(actualV, actualW));

                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], Wtilde[i][j]) + Istim);

                            // Update W explicitly (Heun)
                            W[i][j] = actualW + (deltatODE * reactionW(Vtilde[i][j], Wtilde[i][j]));
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    /*-------------------------
    --  ADI1.5 (1.5 order)  --
    --------------------------*/
    else if (strcmp(method, "ADI1.5") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));

                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], actualW) + Istim);

                            // Update W explicitly
                            W[i][j] = actualW + (deltatODE * reactionW(Vtilde[i][j], actualW));
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    /*-------------------------
    --  ADI1.5 (1.5.2 order)  --
    --------------------------*/
    else if (strcmp(method, "ADI1.5.2") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));

                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], actualW) + Istim);

                            // Update W explicitly
                            W[i][j] = actualW + (deltatODE * reactionW(actualV, actualW));
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    /*-------------------------
    --  ADI1.5 (1.5.3 order)  --
    --------------------------*/
    else if (strcmp(method, "ADI1.5.3") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW, Rw) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));
                            
                            Rw = deltatODE * reactionW(actualV, actualW);
                            Wtilde[i][j] = actualW + 0.5 * Rw;
                            
                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], Wtilde[i][j]) + Istim);

                            // Update W explicitly
                            W[i][j] = actualW + Rw;
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    /*--------------------
    --  Forward Euler   --
    ----------------------*/
    else if (strcmp(method, "FE") == 0)
    {
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];

                            // Update V and W
                            Vtilde[i][j] = actualV + deltatODE * (reactionV(actualV, actualW) + Istim);
                            V[i][j] = Vtilde[i][j];
                            W[i][j] = actualW + deltatODE * reactionW(actualV, actualW);
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        V[i][j] = Vtilde[i][j] + phi * (iDiffusion2nd(i, j, N, Vtilde, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, Vtilde, discFibxMax, discFibxMin, discFibyMax, discFibyMin));
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    // Write infos to file
    char infosFileName[MAX_STRING_SIZE], aux[MAX_STRING_SIZE];
    sprintf(infosFileName, "infos-%d-%.3lf-%.3lf.txt", numberThreads, deltatODE, deltatPDE);
    FILE *fpInfos;
    sprintf(aux, "%s/%s", pathToSaveData, infosFileName);
    fpInfos = fopen(aux, "a");
    fprintf(fpInfos, "Real total execution time: %lf seconds\n", elapsedTotal); 
    

    // Close files
    fclose(fpInfos);

    // Free memory
    free(time);
    free(V);
    free(W);
    free(Vtilde);
    free(Wtilde);
    free(Rv);
    free(rightside);
    free(solution);
    free(c_);
    free(d_);
}

void runMethodConvergenceAnalysis(bool fibro, char *method, float delta_x, float delta_t, int numberThreads, int spatialRate, int rate)
{
    // Update delta x and delta y
    deltax = delta_x;
    deltay = delta_x;

    // Update delta t for ODEs and PDEs
    float deltatODE = delta_t;
    float deltatPDE = delta_t;

    // Get fibro
    bool haveFibrosis = fibro;

    // Number of steps
    int N = round(L / deltax) + 1;                          // Spatial steps (square tissue)
    int M = round(T / deltatPDE) + 1;                       // Number of time steps
    int PdeOdeRatio = round(deltatPDE / deltatODE);     // Ratio between PDE and ODE time steps

    // Allocate and populate time array
    double *time;
    // allocateAndPopulateTime(M, time, deltatPDE);
    time = (double *)malloc(M * sizeof(double));
    for (int i = 0; i < M; i++)
    {
        time[i] = i * deltatPDE;
    }

    // Allocate and initialize variables
    double **V, **W;
    V = (double **)malloc(N * sizeof(double *));
    W = (double **)malloc(N * sizeof(double *));

    for (int i = 0; i < N; i++)
    {
        V[i] = (double *)malloc(N * sizeof(double));
        W[i] = (double *)malloc(N * sizeof(double));
    }
    //initializeVariables(N, V, W);
    initializeVariablesWithSpiral(N, V, W, deltax, rate);

    // Diffusion coefficient - isotropic
    double D = sigma / (chi * Cm);
    double phi = D * deltatPDE / (deltax * deltax);    // For Thomas algorithm - isotropic

    // Variables
    int i, j;                                          // i for y-axis and j for x-axis
    double actualV, actualW;
    double Rw;
    double **Vtilde, **Wtilde, **Rv, **rightside, **solution;
    Vtilde = (double **)malloc(N * sizeof(double *));
    Wtilde = (double **)malloc(N * sizeof(double *));
    Rv = (double **)malloc(N * sizeof(double *));
    rightside = (double **)malloc(N * sizeof(double *));
    solution = (double **)malloc(N * sizeof(double *));

    // Auxiliary arrays for Thomas algorithm 2nd order approximation
    double **c_ = (double **)malloc((N) * sizeof(double *));
    double **d_ = (double **)malloc((N) * sizeof(double *));
    for (int i = 0; i < N; i++)
    {
        Vtilde[i] = (double *)malloc(N * sizeof(double));
        Wtilde[i] = (double *)malloc(N * sizeof(double));
        Rv[i] = (double *)malloc(N * sizeof(double));
        rightside[i] = (double *)malloc(N * sizeof(double));
        solution[i] = (double *)malloc(N * sizeof(double));
        c_[i] = (double *)malloc(N * sizeof(double));
        d_[i] = (double *)malloc(N * sizeof(double));
    }

    // Discretized limits of stimulation area
    int discS1xLimit = round(stim1xLimit / deltax);
    int discS1yLimit = round(stim1yLimit / deltay);
    int discS2xMax = round(stim2xMax / deltax);
    int discS2xMin = round(stim2xMin / deltax);
    int discS2yMax = N;
    int discS2yMin = N - round(stim2yMax / deltay);

    // Discritized limits of fibrotic area
    int discFibxMax = round(fibrosisMaxX / deltax);
    int discFibxMin = round(fibrosisMinX / deltax);
    int discFibyMax = N - round(fibrosisMinY / deltay);
    int discFibyMin = N - round(fibrosisMaxY / deltay);

    // Create directories and files
    char pathToSaveData[MAX_STRING_SIZE];
    if (haveFibrosis)
    {
        createDirectories(pathToSaveData, method, "AFHN-Fibro");
    }
    else
    {
        createDirectories(pathToSaveData, method, "AFHN");
    }


    /*--------------------
    --  ADI 1st order   --
    ----------------------*/
    if (strcmp(method, "ADI1") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V and W without diffusion
                            V[i][j] = actualV + deltatODE * (reactionV(actualV, actualW) + Istim);
                            W[i][j] = actualW + deltatODE * reactionW(actualV, actualW);

                            // Update right side of Thomas algorithm
                            rightside[j][i] = V[i][j];
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines)
                #pragma omp barrier
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, phi, c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, phi, c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        Vtilde[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns)
                #pragma omp barrier
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(Vtilde[i], V[i], N, phi, c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(Vtilde[i], V[i], N, phi, c_[i], d_[i], N, 0);
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }

    
    /*-------------------------
    --  SSI-ADI (2nd order)  --
    --------------------------*/
    else if (strcmp(method, "SSI-ADI") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion and W without diffusion (Heun)
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));
                            Wtilde[i][j] = actualW + (0.5 * deltatODE * reactionW(actualV, actualW));

                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], Wtilde[i][j]) + Istim);

                            // Update W explicitly
                            W[i][j] = actualW + (deltatODE * reactionW(Vtilde[i][j], Wtilde[i][j]));
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    /*-------------------------
    --  ADI1.5 (1.5 order)  --
    --------------------------*/
    else if (strcmp(method, "ADI1.5") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));
                            
                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], actualW) + Istim);

                            // Update W explicitly
                            W[i][j] = actualW + deltatODE * reactionW(Vtilde[i][j], actualW);
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    /*-------------------------
    --  ADI1.5 (1.5.2 order)  --
    --------------------------*/
    else if (strcmp(method, "ADI1.5.2") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));
                            
                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], actualW) + Istim);

                            // Update W explicitly
                            W[i][j] = actualW + deltatODE * reactionW(actualV, actualW);
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    /*-------------------------
    --  ADI1.5 (1.5.3 order)  --
    --------------------------*/
    else if (strcmp(method, "ADI1.5.3") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW, Rw) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];
                            
                            // Update V with diffusion
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));
                            
                            Rw = deltatODE * reactionW(actualV, actualW);
                            Wtilde[i][j] = actualW + 0.5 * Rw;
                            
                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], Wtilde[i][j]) + Istim);

                            // Update W explicitly
                            W[i][j] = actualW + Rw;
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                // 1st: Implicit y-axis diffusion (lines): right side with explicit x-axis diffusion (columns)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibyMin && i <= discFibyMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], discFibxMin, discFibxMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }

                    // Update V
                    for (j = 0; j < N; j++)
                    {
                        V[j][i] = solution[i][j];
                    }
                }

                // 2nd: Implicit x-axis diffusion (columns): right side with explicit y-axis diffusion (lines)
                // Update right side of Thomas algorithm
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V, discFibxMax, discFibxMin, discFibyMax, discFibyMin))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    // Check if i is in fibrotic region
                    if (i >= discFibxMin && i <= discFibxMax)
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], discFibyMin, discFibyMax);
                    }
                    else
                    {
                        ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i], N, 0);
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }



    /*--------------------
    --  Forward Euler   --
    ----------------------*/
    else if (strcmp(method, "FE") == 0)
    {
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        fibrosisFactor, stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, discFibxMax, discFibxMin, discFibyMax, discFibyMin, \
        c_, d_, Vtilde, Wtilde, Rv, rightside, solution)
        {
            while (timeStepCounter < M)
            {
                // Get time step
                timeStep = time[timeStepCounter];

                // Resolve ODEs
                for (int k = 0; k < PdeOdeRatio; k++)
                {
                    #pragma omp for collapse(2)
                    for (i = 0; i < N; i++)
                    {
                        for (j = 0; j < N; j++)
                        {
                            // Stimulus
                            Istim = stimulus(i, j, timeStep, discS1xLimit, discS1yLimit, discS2xMin, discS2xMax, discS2yMin, discS2yMax);

                            // Get actual V and W
                            actualV = V[i][j];
                            actualW = W[i][j];

                            // Update V and W
                            Vtilde[i][j] = actualV + deltatODE * (reactionV(actualV, actualW) + Istim);
                            V[i][j] = Vtilde[i][j];
                            W[i][j] = actualW + deltatODE * reactionW(actualV, actualW);
                        }
                    }
                }

                // Resolve PDEs (Diffusion)
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < N; j++)
                    {
                        V[i][j] = Vtilde[i][j] + phi * (iDiffusion2nd(i, j, N, Vtilde, discFibxMax, discFibxMin, discFibyMax, discFibyMin) + jDiffusion2nd(i, j, N, Vtilde, discFibxMax, discFibxMin, discFibyMax, discFibyMin));
                    }
                }

                // Update time step counter
                #pragma omp master
                {
                    timeStepCounter++;
                }
                #pragma omp barrier
            }
        }

        // Finish measuring total execution time
        finishTotal = omp_get_wtime();
        elapsedTotal = finishTotal - startTotal;
    }


    // Print total execution time
    printf("Total execution time: %lf seconds\n", elapsedTotal);


    // Write last frame to file
    char lastFileName[MAX_STRING_SIZE], aux[MAX_STRING_SIZE], fullLastV[MAX_STRING_SIZE], fullLastW[MAX_STRING_SIZE];
    sprintf(lastFileName, "last-%d-%.5lf.txt", numberThreads, deltatODE);
    sprintf(fullLastV, "fullLastV-%.5lf.txt", deltatODE);
    sprintf(fullLastW, "fullLastW-%.5lf.txt", deltatODE);
    FILE *fpLast, *fpLastV, *fpLastW;
    sprintf(aux, "%s/%s", pathToSaveData, lastFileName);
    fpLast = fopen(aux, "w");
    sprintf(aux, "%s/%s", pathToSaveData, fullLastV);
    fpLastV = fopen(aux, "w");
    sprintf(aux, "%s/%s", pathToSaveData, fullLastW);
    fpLastW = fopen(aux, "w");
    for (int i = 0; i < N; i+=spatialRate)
    {
        for (int j = 0; j < N; j+=spatialRate)
        {
            fprintf(fpLast, "%lf ", V[i][j]);
        }
        fprintf(fpLast, "\n");
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(fpLastV, "%lf ", V[i][j]);
            fprintf(fpLastW, "%lf ", W[i][j]);
        }
        fprintf(fpLastV, "\n");
        fprintf(fpLastW, "\n");
    }


    // Close files
    fclose(fpLast);
    fclose(fpLastV);
    fclose(fpLastW);

    // Free memory
    free(time);
    free(V);
    free(W);
    free(Vtilde);
    free(Wtilde);
    free(Rv);
    free(rightside);
    free(solution);
    free(c_);
    free(d_);
}
#endif // AFHN

#endif // METHODS_H