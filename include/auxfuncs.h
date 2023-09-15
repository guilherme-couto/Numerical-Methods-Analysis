#ifndef AUXFUNCS_H
#define AUXFUNCS_H

#include "includes.h"


// Adapted for 2nd order approximation
void ThomasAlgorithm2nd(double *d, double *solution, unsigned long N, double phi, double *c_, double *d_, int fibroLower, int fibroUpper)
{   
    // Coefficients
    double a = -phi;            // subdiagonal
    double b = 1 + 2 * phi;     // diagonal (1st and last row)
    double c = - 2 * phi;       // superdiagonal

    // Check if it is in the fibrotic region
    if (0 >= fibroLower && 0 <= fibroUpper)
    {
        a = -phi * fibrosisFactor;
        b = 1 + 2 * phi * fibrosisFactor;
        c = - 2 * phi * fibrosisFactor;
    }
    
    // 1st: update auxiliary arrays
    c_[0] = c / b;
    d_[0] = d[0] / b;

    c = -phi;
    
    for (int i = 1; i <= N - 2; i++)
    {
        // Check if it is in the fibrotic region
        if (i >= fibroLower && i <= fibroUpper)
        {
            a = -phi * fibrosisFactor;
            b = 1 + 2 * phi * fibrosisFactor;
            c = -phi * fibrosisFactor;
        }
        else
        {
            a = -phi;
            b = 1 + 2 * phi;
            c = -phi;
        }

        c_[i] = c / (b - a * c_[i - 1]);
        d_[i] = (d[i] - a * d_[i - 1]) / (b - a * c_[i - 1]);
    }
    
    a = - 2 * phi;
    // Check if it is in the fibrotic region
    if (N - 1 >= fibroLower && N - 1 <= fibroUpper)
    {
        a = - 2 * phi * fibrosisFactor;
        b = 1 + 2 * phi * fibrosisFactor;
        c = -phi * fibrosisFactor;
    }
    else
    {
        a = - 2 * phi;
        b = 1 + 2 * phi;
        c = -phi;
    }
    d_[N - 1] = (d[N - 1] - a * d_[N - 2]) / (b - a * c_[N - 2]);

    a = -phi;

    // 2nd: update solution
    solution[N - 1] = d_[N - 1];
    
    for (int i = N - 2; i >= 0; i--)
    {
        // Check if it is in the fibrotic region
        if (i >= fibroLower && i <= fibroUpper)
        {
            a = -phi * fibrosisFactor;
            b = 1 + 2 * phi * fibrosisFactor;
            c = -phi * fibrosisFactor;
        }
        else
        {
            a = -phi;
            b = 1 + 2 * phi;
            c = -phi;
        }
        solution[i] = d_[i] - c_[i] * solution[i + 1];
    }
}

// Adapted for 2nd order approximation
double iDiffusion2nd(int i, int j, int N, double **V, int discFibxMax, int discFibxMin, int discFibyMax, int discFibyMin)
{
    double result = 0.0;
    if (i == 0)
    {
        result = - 2.0*V[i][j] + 2.0*V[i + 1][j]; 
    }
    else if (i == N - 1)
    {
        result = 2.0*V[i - 1][j] - 2.0*V[i][j]; 
    }
    else
    {
        result = V[i - 1][j] - 2.0*V[i][j] + V[i + 1][j];
    }

    if ((i >= discFibyMin && i <= discFibyMax) && (j >= discFibxMin && j <= discFibxMax))
    {
        result *= fibrosisFactor;
    }
    return result;
}

// Adapted for 2nd order approximation
double jDiffusion2nd(int i, int j, int N, double **V, int discFibxMax, int discFibxMin, int discFibyMax, int discFibyMin)
{
    double result = 0.0;
    if (j == 0)
    {
        result = - 2.0*V[i][j] + 2.0*V[i][j + 1]; 
    }
    else if (j == N - 1)
    {
        result = 2.0*V[i][j - 1] - 2.0*V[i][j]; 
    }
    else
    {
        result = V[i][j - 1] - 2.0*V[i][j] + V[i][j + 1];
    }

    if ((i >= discFibyMin && i <= discFibyMax) && (j >= discFibxMin && j <= discFibxMax))
    {
        result *= fibrosisFactor;
    }
    return result;
}

// Normalize voltage for two-dimensional array
void normalizeVoltage2D(int N, double **V)
{
    double V_max = V[0][0];
    double V_min = V[0][0];

    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < N; k++)
        {
            if (V[i][k] > V_max)
            {
                V_max = V[i][k];
            }
            if (V[i][k] < V_min)
            {
                V_min = V[i][k];
            }
        }
    }

    double V_range = V_max - V_min;

    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < N; k++)
        {
            V[i][k] = (V[i][k] - V_min) / V_range;
        }
    }
}

// Normalize voltage for one-dimensional array
void normalizeVoltage1D(int N, double *V)
{
    double V_max = V[0];
    double V_min = V[0];

    for (int i = 0; i < N; i++)
    {
        if (V[i] > V_max)
        {
            V_max = V[i];
        }
        if (V[i] < V_min)
        {
            V_min = V[i];
        }
    }

    double V_range = V_max - V_min;

    for (int i = 0; i < N; i++)
    {
        V[i] = (V[i] - V_min) / V_range;
    }
}

void createDirectories(char *pathToSaveData, char *method, char *cellModel)
{
    char command[] = "mkdir -p";
    char aux[MAX_STRING_SIZE];
    system("mkdir -p ./simulation-files");
    sprintf(pathToSaveData, "%s/%.3lf", "./simulation-files", deltax);
    sprintf(aux, "%s %s", command, pathToSaveData);
    system(aux);
    sprintf(pathToSaveData, "%s/%s", pathToSaveData, cellModel);
    sprintf(aux, "%s %s", command, pathToSaveData);
    system(aux);
    sprintf(pathToSaveData, "%s/%s", pathToSaveData, method);
    sprintf(aux, "%s %s", command, pathToSaveData);
    system(aux);
}

double stimulus(int i, int j, int timeStep, int discS1xLimit, int discS1yLimit, int discS2xMin, int discS2xMax, int discS2yMin, int discS2yMax)
{
   // Stimulus 1
    if (timeStep >= stim1Begin && timeStep <= stim1Begin + stim1Duration && j <= discS1xLimit)
    {
        return stimStrength;
    }
    // Stimulus 2
    else if (timeStep >= stim2Begin && timeStep <= stim2Begin + stim2Duration && j >= discS2xMin && j <= discS2xMax && i >= discS2yMin && i <= discS2yMax)
    {
        return stimStrength;
    }
    return 0.0;
}

void writeToCheckpointFile(float **matrix, int N, double lastCheckpointTime, double timeStep, double deltatODE, double deltatPDE, double deltax, int spatialRate, char *method, char *cellModel, bool haveFibrosis, double fibrosisFactor)
{
    double presentTime = omp_get_wtime();

    FILE *file = fopen("checkpoint.txt", "w");
    fprintf(file, "CHECKPOINT\n");

    fprintf(file, "\nSimultation T: %.3lf ms\n", T);
    fprintf(file, "Simulation L: %.3lf cm\n", L);
    fprintf(file, "Method: %s\n", method);
    fprintf(file, "Cell model: %s\n", cellModel);
    if (haveFibrosis)
    {
        fprintf(file, "Fibrosis factor: %.2lf\n", fibrosisFactor);
    }
    else
    {
        fprintf(file, "Without fibrosis\n");
    }

    fprintf(file, "\nDelta t (ODE): %.5lf ms\n", deltatODE);
    fprintf(file, "Delta t (PDE): %.5lf ms\n", deltatPDE);
    fprintf(file, "Delta x: %.5lf cm\n", deltax);

    fprintf(file, "\nActual time step: %.5lf ms\n", timeStep);
    double porcentageOfTotalTime = (timeStep / T) * 100;
    fprintf(file, "Porcentage of total time: %.2lf %%\n", porcentageOfTotalTime);
    fprintf(file, "Time since last checkpoint: %.3lf seconds\n", presentTime - lastCheckpointTime);
    
    fprintf(file, "\nMatrix (with spatialRate = %d):\n", spatialRate);
    for (int i = 0; i < N; i+=spatialRate)
    {
        for (int j = 0; j < N; j+=spatialRate)
        {
            fprintf(file, "%lf ", matrix[i][j]);
            // Check if NaN
            if (matrix[i][j] != matrix[i][j])
            {
                printf("NaN found at (%d, %d)\n", i, j);
                fprintf(file, "\nNaN found at (%d, %d)\n", i, j);
                exit(1);
            }
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

// TODO: fix this function. Allocation with function call is not working
// void allocateAndPopulateTime(int M, double *time, double deltatPDE)
// {
//     time = (double *)malloc(M * sizeof(double));
//     for (int i = 0; i < M; i++)
//     {
//         time[i] = i * deltatPDE;
//     }
// }


//############################################
//##                                        ##
//##     Adapted FitzHugh-Nagumo (AFHN)     ##
//##                                        ##
//############################################
#if defined(AFHN)
// Initialize variables
void initializeVariables(int N, double **V, double **W)
{
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < N; k++)
        {
            V[i][k] = V_init;
            W[i][k] = W_init;
        }
    }
}

void initializeVariablesWithSpiral(int N, double **V, double **W, double dxActual, int rate)
{
	double dxRef = 0.00125;
    int NRef = round(L / dxRef) + 1;
	
	int MAX_LINE_LENGTH = 20;
	char lineV[MAX_LINE_LENGTH];
    char lineW[MAX_LINE_LENGTH];
	char *ptr;
	
	int counter = 0;
	
	FILE *fileV, *fileW;
	fileV = fopen("./spiral-conditions/spiralV.txt", "r");
    if(fileV == NULL)
    {
		printf("No file for V initial condition\n");
        return;
	}
	fileW = fopen("./spiral-conditions/spiralW.txt", "r");
    if(fileW == NULL)
    {
		printf("No file for W initial condition\n");
        return;
	}
	
	// Read initial conditions
    int a = 0;
    int b = 0;
    for (int i = 0; i < NRef; i++)
    {
        for (int k = 0; k < NRef; k++)
        {
            fgets(lineV, MAX_LINE_LENGTH, fileV);
            fgets(lineW, MAX_LINE_LENGTH, fileW);

            if (i % rate == 0 && k % rate == 0)
            {
                V[a][b] = strtod(lineV, &ptr);
                W[a][b] = strtod(lineW, &ptr);
                b++;
                if (b == N)
                {
                    a++;
                    b = 0;
                }
            }
        }
    }
    
    fclose(fileV);
    fclose(fileW);
}

// Allocate memory
// TODO: fix this function. Allocation with function call is not working
// void allocateMemory(int N, double **V, double **W)
// {
//     V = (double **)malloc(N * sizeof(double *));
//     W = (double **)malloc(N * sizeof(double *));

//     for (int i = 0; i < N; i++)
//     {
//         V[i] = (double *)malloc(N * sizeof(double));
//         W[i] = (double *)malloc(N * sizeof(double));
//     }
// }
#endif // AFHN



//###########################################
//##                                       ##
//##     Minimal Ventricular Model (MV)    ##
//##                                       ##
//###########################################
#if defined(MV)
// Initialize variables
void initializeVariables(int N, double **U, double **V, double **W, double **S)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            U[i][j] = U_init;
            V[i][j] = V_init;
            W[i][j] = W_init;
            S[i][j] = S_init;
        }
    }
}

// Standard Heaviside function
double H(double x)
{
    if (x > 0.0)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

// Convert the dimensionless voltage variable U to dimensions of mV
double rescaleU(double U)
{
    return 85.7*U - 84.0;
}
#endif // MV



//###########################################
//##                                       ##
//##     ten Tusscher 2006 model (TT2)     ##
//##                                       ##
//###########################################
#if defined(TT2)
#if defined(EPI) || defined(M) || defined(ENDO)
// Initialize variables
void initializeVariables(int N, double **V, double **X_r1, double **X_r2, double **X_s, double **m, double **h, double **j, double **d, double **f, double **f2, double **fCass, double **s, double **r, double **Ca_i, double **Ca_SR, double **Ca_SS, double **R_prime, double **Na_i, double **K_i)
{
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < N; k++)
        {
            V[i][k] = V_init;
            X_r1[i][k] = X_r1_init;
            X_r2[i][k] = X_r2_init;
            X_s[i][k] = X_s_init;
            m[i][k] = m_init;
            h[i][k] = h_init;
            j[i][k] = j_init;
            d[i][k] = d_init;
            f[i][k] = f_init;
            f2[i][k] = f2_init;
            fCass[i][k] = fCass_init;
            s[i][k] = s_init;
            r[i][k] = r_init;
            Ca_i[i][k] = Ca_i_init;
            Ca_SR[i][k] = Ca_SR_init;
            Ca_SS[i][k] = Ca_SS_init;
            R_prime[i][k] = R_prime_init;
            Na_i[i][k] = Na_i_init;
            K_i[i][k] = K_i_init;
        }
    }
}
#endif // EPI || M || ENDO
#endif // TT2

#endif // AUXFUNCS_H