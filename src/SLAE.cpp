#include "SLAE.hpp"

void SLAE::IterativeMethod(int NumOfW)
{
    double *curX, *prevX;
    VecotorCopy(x0, local_x0);

    if (method == 0)
    {
        prevX = local_x0;
        curX = x;
    }
    else
    {
        curX = local_x0;
        prevX = local_x0;
    }

    double normB = VectorNorm(b);
    double RelativeDiscrepancy = CalculateRelativeDiscrepancy(x0);
    int curIteration = 0;
    double DiscrepancyF_Ax = 0;

    for (; curIteration < MaxNumOfIterations && RelativeDiscrepancy > AccuracyOfTheSolution; curIteration++)
    {
        DiscrepancyF_Ax = 0;

        for (int i = 0; i < n; i++)
        {
            int indX = 0;
            double sum = 0;

            for (int j = 0; j < 4; j++)
            {
                indX = i + TableOfShifts[j];
                if (indX >= 0)
                {
                    sum += prevX[indX] * matrix[i][j];
                }
            }
            sum += prevX[i] * matrix[i][4];

            for (int j = 5; j < 9; j++)
            {
                indX = i + TableOfShifts[j];
                if (indX < n)
                {
                    sum += prevX[indX] * matrix[i][j];
                }
            }

            curX[i] = prevX[i] + w * (b[i] - sum) / matrix[i][4];
            DiscrepancyF_Ax += (b[i] - sum) * (b[i] - sum);
        }

        std::swap(curX, prevX);
        RelativeDiscrepancy = sqrt(DiscrepancyF_Ax) / normB;

        if (isinf(RelativeDiscrepancy) || isnan(RelativeDiscrepancy))
            break;
    }

    VecotorCopy(prevX, x);
    TableOfNumOfConditionality[NumOfW - 1] = CalculateNumOfConditionality(RelativeDiscrepancy);

    if (curIteration >= MaxNumOfIterations)
    {
        NumOfIterationsDependingOnW[NumOfW - 1] = -1;
    }
    else if (RelativeDiscrepancy < AccuracyOfTheSolution)
    {
        NumOfIterationsDependingOnW[NumOfW - 1] = curIteration;
        VectorOutput(x);
    }
    else if (isnan(RelativeDiscrepancy) || isinf(RelativeDiscrepancy))
    {
        NumOfIterationsDependingOnW[NumOfW - 1] = -2;
    }
}

void SLAE::Input(FILE *matrixFile, FILE *vectorFile, FILE *paramFile)
{
    fscanf_s(matrixFile, "%d", &n);
    fscanf_s(matrixFile, "%d", &m);

    fscanf_s(paramFile, "%d", &method);
    fscanf_s(paramFile, "%lf", &AccuracyOfTheSolution);
    fscanf_s(paramFile, "%d", &MaxNumOfIterations);

    AllocateMemory();
    InitializeShiftsTable();

    for (int Icount = 0; Icount < 9; Icount++)
    {
        int curDiag = TableOfShifts[Icount];
        if (curDiag <= 0)
        {
            for (int i = abs(curDiag); i < n; i++)
            {
                fscanf_s(matrixFile, "%lf", &matrix[i][Icount]);
            }
        }
        else
        {
            for (int i = 0; i < (n - abs(curDiag)); i++)
            {
                fscanf_s(matrixFile, "%lf", &matrix[i][Icount]);
            }
        }
    }

    for (int i = 0; i < n; i++)
        fscanf_s(vectorFile, "%lf", &b[i]);
    for (int i = 0; i < n; i++)
        fscanf_s(vectorFile, "%lf", &x0[i]);

    for (int i = 0; i < n; i++)
    {
        xtrue[i] = static_cast<double>(i + 1);
    }
}

void SLAE::AllocateMemory()
{
    matrix = new double *[n];
    for (int i = 0; i < n; ++i)
    {
        matrix[i] = new double[n];
    }
    b = new double[n];
    x = new double[n];
    x0 = new double[n];
    xtrue = new double[n];
    local_x0 = new double[n];
    vectorForDiscrepancy = new double[n];
    NumOfIterationsDependingOnW = new int[200];
    TableOfNumOfConditionality = new double[200];
}

void SLAE::MatrixVectorMultiplicationForDiscrepancy(double *first)
{
    for (int i = 0; i < n; i++)
    {
        int indX = 0;
        double sum = 0;

        for (int j = 5; j < 9; j++)
        {
            indX = i + TableOfShifts[j];
            if (indX < n)
            {
                sum += first[indX] * matrix[i][j];
            }
        }
        sum += first[i] * matrix[i][4];

        for (int j = 0; j < 4; j++)
        {
            indX = i + TableOfShifts[j];
            if (indX >= 0)
            {
                sum += first[indX] * matrix[i][j];
            }
        }
        vectorForDiscrepancy[i] = b[i] - sum;
    }
}

double SLAE::CalculateRelativeDiscrepancy(double *first)
{
    MatrixVectorMultiplicationForDiscrepancy(first);
    return VectorNorm(vectorForDiscrepancy) / VectorNorm(b);
}

double SLAE::CalculateNumOfConditionality(double RelativeDiscrepancy)
{
    VecotorSubtract(x, xtrue);
    double VectorXRelDiscrepancy = VectorNorm(vectorForDiscrepancy) / VectorNorm(xtrue);
    return VectorXRelDiscrepancy / RelativeDiscrepancy;
}

void SLAE::OutputDense()
{
    double **matrixDense = new double *[n];
    for (int i = 0; i < n; ++i)
    {
        matrixDense[i] = new double[n];
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            matrixDense[i][j] = 0.0;
        }
    }

    for (int Icount = 0; Icount < 9; Icount++)
    {
        int curDiagonal = TableOfShifts[Icount];
        if (curDiagonal <= 0)
        {
            for (int i = abs(curDiagonal), j = 0; i < n && j < n; i++, j++)
            {
                matrixDense[i][j] = matrix[i][Icount];
            }
        }
        else
        {
            for (int i = 0, j = abs(curDiagonal); i < n && j < n; i++, j++)
            {
                matrixDense[i][j] = matrix[i][Icount];
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf(REALOUTD, matrixDense[i][j]);
        }
        printf("\n");
    }
}

double SLAE::VectorNorm(double *first)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
    {
        norm += first[i] * first[i];
    }
    return sqrt(norm);
}

void SLAE::VecotorSubtract(double *first, double *second)
{
    for (int i = 0; i < n; i++)
    {
        vectorForDiscrepancy[i] = first[i] - second[i];
    }
}

void SLAE::VecotorCopy(double *first, double *second)
{
    for (int i = 0; i < n; i++)
    {
        second[i] = first[i];
    }
}

void SLAE::InitializeShiftsTable()
{
    for (int i = 0; i < 3; i++)
    {
        TableOfShifts[i] -= m;
        TableOfShifts[8 - i] += m;
    }
}

void SLAE::OutputSolutionVector(FILE *out)
{
    for (int i = 0; i < n; i++)
    {
        fprintf_s(out, REALOUT, x[i]);
    }
    fprintf_s(out, "\n");
}

void SLAE::VectorOutput(double *curX)
{
    for (int i = 0; i < n; i++)
    {
        printf_s(REALOUT, curX[i]);
    }
    printf_s("\n");
}

void SLAE::OutputResultParametrs()
{
    for (int i = 1; i <= 200; i++)
    {
        printf("%.2lf ", 0.01 * i);
        printf("%d ", NumOfIterationsDependingOnW[i - 1]);
        printf("%lf\n", TableOfNumOfConditionality[i - 1]);
    }
}
