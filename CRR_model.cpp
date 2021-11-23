#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <windows.h>

LARGE_INTEGER timeStart;
LARGE_INTEGER timeEnd;
LARGE_INTEGER frequency;

clock_t start, end;
double payoff[10000];
double tused;

int main()
{
    QueryPerformanceFrequency(&frequency);
    double quadpart = (double)frequency.QuadPart;

    double S = 5, X = 5, T = 1, r = 0.15, Vol = 0.25, q = 0;
    double u, d, Pu, Pd, dT, current;
    FILE *Write = fopen("CRR_Result.txt", "w");

    int n = 10000;
    dT = (double)(T / n);
    u = exp(Vol * sqrt(dT));
    d = exp(-Vol * sqrt(dT));
    Pu = (exp((r - q) * dT) - d) / (u - d);
    Pd = 1 - Pu;
    printf("%lf, %lf, %lf, %lf, %lf\n", dT, u, d, Pu, Pd);

    //計算最後一期Payoff
    current = S * pow(u, n);
    for (int i = 0; i <= n; i++)
    {
        // payoff[i] = fmax(current - X, 0.0);
        payoff[i] = fmax((current - 4) * (current - 5) * (current - 6) * (current - 7) + 5 - X, 0.0);
        current *= d * d;
    }
    //Backward induction
    for (int i = n; i > 0; i--)
    {
        current = S * pow(u, i - 1);
        for (int j = 0; j < i; j++)
        {
            // payoff[j] = fmax(current - X, exp(-r * dT) * (Pu * payoff[j] + Pd * payoff[j + 1]));
            payoff[j] = fmax((current - 4) * (current - 5) * (current - 6) * (current - 7) + 5 - X, exp(-r * dT) * (Pu * payoff[j] + Pd * payoff[j + 1]));
            current *= d * d;
        }
    }
    double mean = payoff[0];
    // printf("%lf", mean);

    for (int n = 200; n <= 5000; n += 10)
    {
        QueryPerformanceCounter(&timeStart);
        //計算CRR相關參數
        dT = (double)T / n;
        u = exp(Vol * sqrt(dT));
        d = exp(-Vol * sqrt(dT));
        Pu = (exp((r - q) * dT) - d) / (u - d);
        Pd = 1 - Pu;
        //計算最後一期Payoff
        current = S * pow(u, n);
        for (int i = 0; i <= n; i++)
        {
            payoff[i] = fmax((current - 4) * (current - 5) * (current - 6) * (current - 7) + 5 - X, 0);
            current *= d * d;
        }
        //Backward induction
        for (int i = n; i > 0; i--)
        {
            current = S * pow(u, i - 1);
            for (int j = 0; j < i; j++)
            {
                payoff[j] = fmax((current - 4) * (current - 5) * (current - 6) * (current - 7) + 5 - X, exp(-r * dT) * (Pu * payoff[j] + Pd * payoff[j + 1]));
                current *= d * d;
            }
        }
        QueryPerformanceCounter(&timeEnd);
        tused = (timeEnd.QuadPart - timeStart.QuadPart) / quadpart;
        fprintf(Write, "%d, %lf,%lf, %lf\n", n, payoff[0], fabs(payoff[0] - mean), tused);
    }
}