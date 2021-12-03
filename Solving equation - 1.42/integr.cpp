//
// Created by osboxes on 12/2/21.
//
//Code from https://github.com/raxtemur/archive/blob/master/Projects/NM%20Integration/main.c

#include "integr.h"
#include <math.h>
#include <stdio.h>

#define EPS 0.000000001



int countSteps(integral I, function fx, double a, double b, int m, double eps, unsigned long long *ans)
{
    unsigned long long Im;
    double I1, I0;
    double p;
    double eps0, h = (b - a)/32;
    if(eps < EPS)
    {
        printf("internal error\n");
        return(0);
    }
    if(!Integral_midPoint(fx, a, b, eps, 16, &I0))
    {
        printf("internal error\n");
        return(0);
    }
    if(!Integral_midPoint(fx, a, b, eps, 32, &I1))
    {
        printf("internal error\n");
        return(0);
    }
    eps0 = (I1 - I0)/(pow(2, m + 1) - 1);
    if (eps0 > eps)
    {
        h = h*pow(eps/eps0, 1.0/(m + 1));
        *ans = ceil(fabs((b - a) / h));
    }
    else
    {
        *ans = 32;
    }
    return 1;
}

//method works if it's known that 2nd derivative is bounded(better <1)
int countSteps_midPoint(function fx, double a, double b, int m, double eps, unsigned long long *ans)
{
    *ans = round(sqrt(pow(abs(b-a), (m+1))/eps)) + 1;
    return 1;
}

int Integral_midPoint(function fx, double a, double b, double eps, unsigned long long steps, double *ans)//m = 1, p = 2;
{
    double h, sum = 0;
    double x, f;
    if (steps == 0)
    {
        countSteps_midPoint(fx, a, b, 1, eps, &steps);
        /*if(!countSteps_midPoint(fx, a, b, 1, eps, &steps))
        {
            printf("internal error\n");
            return(0);
        }*/

    }
    printf("steps: %ld\n", steps);
    h = (b - a)/steps;
    x = a;
    for(int i = 0; i < steps; i++)
    {
        fx((x + h/2), &f);
        sum += f;
        x = x + h;
    }
    *ans = sum*h;
    return 1;
}

int Integral_Gauss(function fx, double a, double b, double eps, unsigned long long steps, double *ans)//m = 2
{
    double h, sum = 0;
    double x, x1, x2, f1, f2;
    if (steps == 0)
    {
        if(!countSteps(Integral_Gauss, fx, a, b, 2, eps, &steps))
        {
            printf("internal error\n");
            return(0);
        }

    }
    h = (b - a)/steps;
    x = a;
    for(int i = 0; i < steps; i++)
    {
        x1 = (x + h/2) - (h/(2*sqrt(3)));
        x2 = (x + h/2) + (h/(2*sqrt(3)));
        fx(x1, &f1);
        fx(x2, &f2);
        sum += (f1 + f2)*(h/2);
        x += h;

    }
    *ans = sum;
    return 1;
}

int Integral_Simpson(function fx, double a, double b, double eps, unsigned long long steps, double *ans)//m = 2
{
    double h, sum = 0;
    double x, f1, f2, f3;
    if (steps == 0)
    {
        if(!countSteps(Integral_Simpson, fx, a, b, 2, eps, &steps))
        {
            printf("internal error\n");
            return(0);
        }

    }
    h = (b - a)/steps;
    x = a;
    for(int i = 0; i < steps; i++)
    {
        fx(x, &f1);
        fx(x + h/2, &f2);
        fx(x + h, &f3);
        sum += h*(f1 + 4*f2 + f3)/6;
        x = x + h;
    }
    *ans = sum;
    return 1;
}

//int countSteps(integral I, function fx, double a, double b, int m, double eps, unsigned long long *ans);
//int Integral_midPoint(function fx, double a, double b, double eps, unsigned long long steps, double *ans);
