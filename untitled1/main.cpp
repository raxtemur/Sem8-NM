#include <iostream>
#include <stdio.h>
#include <math.h>
#include "integr.h"

using namespace std;

int f(double x, double *ans)
{
    *ans = sqrt(x*x*x+1/x);
    return 1;
}

int fun(double x, double alpha, double norm_coeff, double *ans)
{
    double a = 0, b = x, eps = 0.0001;
    //Integral_Gauss(f, a, b, eps, 0, ans);
    Integral_midPoint(f, a, b, eps, 0, ans);
    //Integral_Simpson(f, a, b, eps, 0, ans);
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
    while (abs((x0 - x1)/x0) > 1e-5)
    {
        x0 = x1;
        fun(x0, alpha, norm_coeff, &x1);
        x1 = x0 - x1;

        printf("%lf : %lf\n", x0, x1);
    }
    *ans = x1;
}

int main() {
    double a, b, eps = 0.0001, ans, alpha , norm_c;
    cout << "Input alpha:";
    cin >> alpha;
    int NumOfSign = ceil(-log(eps)/log(10) + 1);

    fun(10, alpha, 1, &ans);
    printf("f(10)=%lf\n", ans);


    segment(alpha, &a, &b, &norm_c);
    printf("limits: %lf, %lf\n", a, b);
    iterate(a, b, alpha, norm_c, &ans);
    //fun(10, alpha, 1, &ans);
    printf("%.*lf\n", NumOfSign, ans);

    return 0;
}
