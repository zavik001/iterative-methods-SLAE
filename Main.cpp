#include <iostream>
#include <math.h>
#include "SLAE.hpp"

int main()
{
    FILE *inputMatrix, *inputVector, *inputParam, *outFile;

    fopen_s(&inputMatrix, "data/matrix.txt", "r");
    fopen_s(&inputVector, "data/vector.txt", "r");
    fopen_s(&inputParam, "data/param.txt", "r");
    fopen_s(&outFile, "data/out.txt", "w+");

    SLAE slae;
    slae.Input(inputMatrix, inputVector, inputParam);
    slae.OutputDense();

    if (slae.method == 0) // Якоби
    {
        for (int i = 1; i <= 200; i++)
        {
            slae.w += 0.01;
            slae.IterativeMethod(i);
            slae.OutputSolutionVector(outFile);
        }
    }
    else // Зейдель
    {
        for (int i = 1; i <= 200; i++)
        {
            slae.w += 0.01;
            slae.IterativeMethod(i);
            slae.OutputSolutionVector(outFile);
        }
    }
    slae.OutputResultParametrs();
    return 0;
}
