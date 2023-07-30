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

// For velocity
bool S1VelocityTag = true;
double S1Velocity = 0.0;


//############################################
//##                                        ##
//##     Adapted FitzHugh-Nagumo (AFHN)     ##
//##                                        ##
//############################################
#if defined(AFHN)
void runMethod(char *options[], float deltatODE, float deltatPDE, int numberThreads)
{
    // Get options
    char *method = options[0];
    char *measureTotalTime = options[1];
    char *saveDataToError = options[2];
    char *saveDataToGif = options[3];
    char *measureTimeParts = options[4];
    char *measureS1Velocity = options[5];
    char *evaluateVulnerabilityWindow = options[6];

    // Number of steps
    int N = (int)(L / deltax);                          // Spatial steps (square tissue)
    int M = (int)(T / deltatPDE);                       // Number of time steps
    int PdeOdeRatio = (int)(deltatPDE / deltatODE);     // Ratio between PDE and ODE time steps

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
    int i, j;
    double actualV, actualW;
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
    int discS1xLimit = (int)(stim1xLimit / deltax);
    int discS1yLimit = (int)(stim1yLimit / deltay);
    int discS2xMax = (int)(stim2xMax / deltax);
    int discS2xMin = (int)(stim2xMin / deltax);
    int discS2yMax = N;
    int discS2yMin = N - (int)(stim2yMax / deltay);

    // File names
    char framesFileName[MAX_STRING_SIZE], lastFrameFileName[MAX_STRING_SIZE], infosFileName[MAX_STRING_SIZE];
    sprintf(framesFileName, "frames-%d-%.3lf-%.3lf.txt", numberThreads, deltatODE, deltatPDE);
    sprintf(lastFrameFileName, "last-%d-%.3lf-%.3lf.txt", numberThreads, deltatODE, deltatPDE);
    sprintf(infosFileName, "infos-%d-%.3lf-%.3lf.txt", numberThreads, deltatODE, deltatPDE);
    int saverate = ceil(M / 100.0);
    FILE *fpFrames, *fpLast, *fpInfos;

    // Create directories and files
    char pathToSaveData[MAX_STRING_SIZE];
    createDirectories(pathToSaveData, method, "AFHN");
    
    // File pointers
    char aux[MAX_STRING_SIZE];
    sprintf(aux, "%s/%s", pathToSaveData, framesFileName);
    fpFrames = fopen(aux, "w");
    sprintf(aux, "%s/%s", pathToSaveData, lastFrameFileName);
    fpLast = fopen(aux, "w");
    sprintf(aux, "%s/%s", pathToSaveData, infosFileName);
    fpInfos = fopen(aux, "w");

    /*--------------------
    --  ADI 1st order   --
    ----------------------*/
    if (strcmp(method, "ADI1") == 0)
    {   
        // Start measuring total execution time
        startTotal = omp_get_wtime();

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, \
        c_, d_, Vtilde, Wtilde, S1Velocity, S1VelocityTag, saverate, fpFrames, fpLast, fpInfos, \
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
                    ThomasAlgorithm2nd(rightside[i], solution[i], N, phi, c_[i], d_[i]);

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
                    ThomasAlgorithm2nd(Vtilde[i], V[i], N, phi, c_[i], d_[i]);
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
                    // Write frames to file
                    if (timeStepCounter % saverate == 0)
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
                            fprintf(fpInfos, "S1 velocity: %lf m/s\n", S1Velocity);
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

        #pragma omp parallel num_threads(numberThreads) default(none) private(i, j, Istim, actualV, actualW) \
        shared(V, W, N, M, L, T, D, phi, deltatODE, deltatPDE, time, timeStep, timeStepCounter, PdeOdeRatio, \
        stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, \
        c_, d_, Vtilde, Wtilde, S1Velocity, S1VelocityTag, saverate, fpFrames, fpLast, fpInfos, \
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
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V) + jDiffusion2nd(i, j, N, V))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));
                            Wtilde[i][j] = actualW + (0.5 * deltatODE * reactionW(actualV, actualW));

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
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i]);

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
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i]);
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
                    // Write frames to file
                    if (timeStepCounter % saverate == 0)
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
                            fprintf(fpInfos, "S1 velocity: %lf m/s\n", S1Velocity);
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
        stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, \
        c_, d_, Vtilde, Wtilde, S1Velocity, S1VelocityTag, saverate, fpFrames, fpLast, fpInfos, \
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
                            Vtilde[i][j] = actualV + ((0.5 * phi) * (iDiffusion2nd(i, j, N, V) + jDiffusion2nd(i, j, N, V))) + (0.5 * deltatODE * (reactionV(actualV, actualW) + Istim));

                            // Update V reaction term
                            Rv[i][j] = 0.5 * deltatODE * (reactionV(Vtilde[i][j], actualW) + Istim);

                            // Update W explicitly (Heun)
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
                        rightside[j][i] = (V[i][j] + (0.5 * phi * jDiffusion2nd(i, j, N, V))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    ThomasAlgorithm2nd(rightside[i], solution[i], N, (0.5 * phi), c_[i], d_[i]);

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
                        rightside[i][j] = (V[i][j] + (0.5 * phi * iDiffusion2nd(i, j, N, V))) + Rv[i][j];
                    }
                }

                // Solve tridiagonal system for V (Thomas algorithm)
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    ThomasAlgorithm2nd(rightside[i], V[i], N, (0.5 * phi), c_[i], d_[i]);
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
                    // Write frames to file
                    if (timeStepCounter % saverate == 0)
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
                            fprintf(fpInfos, "S1 velocity: %lf m/s\n", S1Velocity);
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
        stimStrength, stim1Duration, stim2Duration, stim1Begin, stim2Begin, stim1xLimit, stim1yLimit, \
        discS1xLimit, discS1yLimit, discS2xMax, discS2yMax, discS2xMin, discS2yMin, \
        c_, d_, Vtilde, Wtilde, S1Velocity, S1VelocityTag, saverate, fpFrames, fpLast, fpInfos, \
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
                    for (i = 1; i < N - 1; i++)
                    {
                        for (j = 1; j < N - 1; j++)
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

                // Boundary conditions
                #pragma omp for
                for (i = 0; i < N; i++)
                {
                    Vtilde[i][0] = Vtilde[i][1];
                    Vtilde[i][N-1] = Vtilde[i][N-2];
                    Vtilde[0][i] = Vtilde[1][i];
                    Vtilde[N-1][i] = Vtilde[N-2][i];
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
                for (i = 1; i < N - 1; i++)
                {
                    for (j = 1; j < N - 1; j++)
                    {
                        // Update V
                        V[i][j] = Vtilde[i][j] + phi * (Vtilde[i-1][j] - 2.0*Vtilde[i][j] + Vtilde[i+1][j]);
                        V[i][j] += phi * (Vtilde[i][j-1] - 2.0*Vtilde[i][j] + Vtilde[i][j+1]);
                    }
                }

                // Boundary conditions
                #pragma omp for
                for (i = 0; i < N; i++)
                {
                    V[i][0] = V[i][1];
                    V[i][N-1] = V[i][N-2];
                    V[0][i] = V[1][i];
                    V[N-1][i] = V[N-2][i];
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
                    // Write frames to file
                    if (timeStepCounter % saverate == 0)
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
                            fprintf(fpInfos, "S1 velocity: %lf m/s\n", S1Velocity);
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
    fprintf(fpInfos, "Total execution time: %lf seconds\n", elapsedTotal);
    fprintf(fpInfos, "PDE/ODE ratio: %d\n", PdeOdeRatio);
    fprintf(fpInfos, "ODE execution time: %lf seconds\n", elapsedODE);
    fprintf(fpInfos, "PDE execution time: %lf seconds\n", elapsedPDE);

    // Write last frame to file
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(fpLast, "%lf ", V[i][j]);
        }
        fprintf(fpLast, "\n");
    }

    // Close files
    fclose(fpFrames);
    fclose(fpLast);
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
#endif // AFHN

#endif // METHODS_H