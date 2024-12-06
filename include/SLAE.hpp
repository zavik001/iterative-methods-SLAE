#include <cstdio>
#include <math.h>
#include <iostream>

#define REALOUT "%.15lf\n"
#define REALOUTD "%.0lf\t"
#define EPS 1e-13

using namespace std;

class SLAE
{
public:
    int n, m, MaxNumOfIterations, method;
    int *NumOfIterationsDependingOnW;
    double AccuracyOfTheSolution, w = 0.00;
    double *x, *x0, *b, *vectorForDiscrepancy, *local_x0, *TableOfNumOfConditionality, *xtrue;
    double **matrix;
    int TableOfShifts[9] = {-4, -3, -2, -1, 0, 1, 2, 3, 4};

    void Input(FILE *matrixFile, FILE *vectorFile, FILE *paramFile);

    void OutputDense();

    void IterativeMethod(int NumOfW);

    void VectorOutput(double *curX);

    void OutputSolutionVector(FILE *out);

    void InitializeShiftsTable();

    double VectorNorm(double *first);

    double CalculateRelativeDiscrepancy(double *first);

    double CalculateNumOfConditionality(double RelativeDiscrepancy);

    void MatrixVectorMultiplicationForDiscrepancy(double *vectorMult);

    void VecotorCopy(double *first, double *second);

    void VecotorSubtract(double *first, double *second);

    void OutputResultParametrs();

protected:
    void AllocateMemory();
};