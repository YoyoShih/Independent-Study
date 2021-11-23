#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <windows.h>

LARGE_INTEGER timeStart;
LARGE_INTEGER timeEnd;
LARGE_INTEGER frequency;

double tused;
clock_t start, end;
double Matrix[200000];

int main()
{
    QueryPerformanceFrequency(&frequency);
    double quadpart = (double)frequency.QuadPart;

    double S = 5, X = 5, T = 1, r = 0.15, Sigma = 0.25, q = 0, R;

    double u, d, Pu, Pd, dT, current;
    int n = 10000;
    dT = (double)T / n;
    u = exp(Sigma * sqrt(dT));
    d = exp(-Sigma * sqrt(dT));
    Pu = (exp((r - q) * dT) - d) / (u - d);
    Pd = 1.0 - Pu;
    //計算最後一期Payoff
    current = S * pow(u, n);
    for (int i = 0; i <= n; i++)
    {
        // Matrix[i] = fmax(current - X, 0);
        Matrix[i] = fmax((current - 4) * (current - 5) * (current - 6) * (current - 7) + 5 - X, 0.0);
        // Matrix[i] = fmax((current - 0.7 * X) * (current - 1.3 * S) - 5000, 0);
        current *= d * d;
    }
    //Backward induction
    for (int i = n; i > 0; i--)
    {
        current = S * pow(u, i - 1);
        for (int j = 0; j < i; j++)
        {
            // Matrix[j] = fmax(current - X, exp(-r * dT) * (Pu * Matrix[j] + Pd * Matrix[j + 1]));
            // Matrix[j] = fmax((current - 0.7 * X) * (current - 1.3 * S) - 5000, exp(-r * dT) * (Pu * Matrix[j] + Pd * Matrix[j + 1]));
            Matrix[j] = fmax((current - 4) * (current - 5) * (current - 6) * (current - 7) + 5 - X, exp(-r * dT) * (Pu * Matrix[j] + Pd * Matrix[j + 1]));
            current *= d * d;
        }
    }
    double mean = Matrix[0];
    // printf("%lf", mean);

    double Pm, DeltaT;
    FILE *Write = fopen("KRL_Result.txt", "w");
    for (int n = 200; n <= 5000; n += 10)
    {
        QueryPerformanceCounter(&timeStart);
        //count start
        DeltaT = (double)T / n;
        double lambda = 1.224745;
        R = exp(r * DeltaT);
        u = exp(lambda * Sigma * sqrt(DeltaT));
        d = exp(-lambda * Sigma * sqrt(DeltaT));
        Pu = 1 / (2 * lambda * lambda) + (r - q - Sigma * Sigma / 2.0) * sqrt(DeltaT) / (2 * lambda * Sigma);
        Pd = 1 / (2 * lambda * lambda) - (r - q - Sigma * Sigma / 2.0) * sqrt(DeltaT) / (2 * lambda * Sigma);
        Pm = 1.0 - Pu - Pd;

        //American call option
        for (int i = 0; i <= 2 * n; i++)
        {
            double St = S * pow(u, n - i);
            // Matrix[i] = fmax(0, St - X);
            // Matrix[i] = fmax(0, (St - 0.7 * X) * (St - 1.3 * S) - 5000);
            Matrix[i] = fmax((St - 4) * (St - 5) * (St - 6) * (St - 7) + 5 - X, 0);
        }
        for (int j = n - 1; j >= 0; j--)
        {
            for (int i = 0; i <= 2 * j; i++)
            {
                double St = S * pow(u, j - i);
                Matrix[i] = (Pu * Matrix[i] + Pm * Matrix[i + 1] + Pd * Matrix[i + 2]) / R;
                // Matrix[i] = fmax(Matrix[i], St - X);
                // Matrix[i] = fmax(Matrix[i], (St - 0.7 * X) * (St - 1.3 * S) - 5000);
                Matrix[i] = fmax(Matrix[i], (St - 4) * (St - 5) * (St - 6) * (St - 7) + 5 - X);
            }
        }
        // printf("%lf", Matrix[0]);
        QueryPerformanceCounter(&timeEnd);
        tused = (timeEnd.QuadPart - timeStart.QuadPart) / quadpart;
        // fprintf(Write, "%d, %lf,%lf,%lf, %lf\n", n, Matrix[0], Matrix[0] - mean, fabs(Matrix[0] - mean), tused);
        fprintf(Write, "%d, %lf,%lf,%lf\n", n, Matrix[0], fabs(Matrix[0] - mean), tused);
    }
    // printf("%lf", Matrix[0]);
    return 0;
}
