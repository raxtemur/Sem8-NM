#include <iostream>
#include <stdio.h>
#include <math.h>
#include "integr.h"

#define eps  0.0000001

using namespace std;

int f(double x, double *ans)
{
    *ans = sqrt(x*x*x+1/x);
    return 1;
}

int fun(double x, double alpha, double norm_coeff, double *ans)
{
    double a = 0.1, b = x;
    const double I1 = 0.63245904562901;
    //Integral_Gauss(f, a, b, eps, 0, ans);
    Integral_midPoint(f, a, b, eps, 0, ans);
    //Integral_Simpson(f, a, b, eps, 0, ans);
    *ans += I1;
    *ans -= (5*sqrt(x) + alpha*alpha);
    *ans *= norm_coeff;
}

int segment(double alpha, double *a, double *b, double *norm_coeff)
{
    double ans, step = 0.1, x = step;
    fun(x, alpha, 1, &ans);
    while (ans < 0) {
        x += step;
        fun(x, alpha, 1, &ans);
    } ;
    *a = x - step;
    *b = x + step;
    *norm_coeff = 1/(sqrt((*b)*(*b)*(*b))+1);
    return 0;
}

int iterate(double a, double b, double alpha, double norm_coeff, double *ans)
{
    double x0 = a, x1 = (b+a)/2;
    while (abs((x0 - x1)/x0) > eps)
    {
        x0 = x1;
        fun(x0, alpha, norm_coeff, &x1);
        x1 = x0 - x1;

        printf("%lf : %lf\n", x0, x1);
    }
    *ans = x1;
}

int main() {
    double a, b, ans, alpha , norm_c;
    cout << "Input alpha:";
    cin >> alpha;
    int NumOfSign = ceil(-log(eps)/log(10) + 1);

    //Integral_midPoint(f, 0.1, 10, eps, 0, &ans);
    //ans+=0.63245904562901;
    //printf("\ng(10)=%.*lf\n\n", NumOfSign, ans);

    //fun(10, alpha, 1, &ans);
    //printf("\nf(10)=%.*lf\n\n", NumOfSign, ans);


    segment(alpha, &a, &b, &norm_c);
    printf("limits: %lf, %lf\n", a, b);
    iterate(a, b, alpha, norm_c, &ans);
    //fun(10, alpha, 1, &ans);
    printf("%.*lf\n", NumOfSign, ans);

    return 0;
}
