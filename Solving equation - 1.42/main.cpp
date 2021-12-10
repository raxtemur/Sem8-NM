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
    Integral_midPoint(f, a, b, 0.1*eps, 0, ans);
    //Integral_Simpson(f, a, b, eps, 0, ans);
    *ans += I1;
    *ans -= (5*sqrt(x) + alpha*alpha);
    *ans *= norm_coeff;
}

int fun_p(double x, double *ans)
{
    *ans = sqrt(x*x*x+1/x) - 2.5/sqrt(x);
    return 1;
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
    fun_p(*b, norm_coeff);
    *norm_coeff = 1/(*norm_coeff);//1/(sqrt((*b)*(*b)*(*b))+1);
    return 0;
}

//note, the derivative of fun is negative! So if current x is greater than the root, the next one will be smaller.
//That is the reason why it is always |x0 - x1| > |x' - x1|, where x' is real root. So |x0 - x1| could be used in stopping criteria.
int iterate(double a, double b, double alpha, double norm_coeff, double *ans, FILE *fp)
{
    double x0 = a, x1 = (b+a)/2;
    double q;
    int N;

    fun_p(a, &q);
    q = 1 - q*norm_coeff;
    N = round(log(a*10e-5/(b-a))/log(q))+1;

    for (int i = 0; i<N; i++)
    {
        x0 = x1;
        fun(x0, alpha, norm_coeff, &x1);
        x1 = x0 - x1;

        fprintf(fp, "%lf : %lf\n", x0, x1);
    }
    *ans = x1;
}

int main() {
    double a, b, ans, alpha , norm_c;
    cout << "Input alpha:";
    cin >> alpha;
    int NumOfSign = ceil(-log(eps)/log(10) + 1);
    char filename[32];
    snprintf(filename, sizeof(filename), "Output_%lf.txt", alpha);
    FILE *fp = fopen(filename, "w");
    //Integral_midPoint(f, 0.1, 10, eps, 0, &ans);
    //ans+=0.63245904562901;
    //printf("\ng(10)=%.*lf\n\n", NumOfSign, ans);

    //fun(10, alpha, 1, &ans);
    //printf("\nf(10)=%.*lf\n\n", NumOfSign, ans);


    segment(alpha, &a, &b, &norm_c);
    fprintf(fp,"limits: %lf, %lf\n", a, b);
    iterate(a, b, alpha, norm_c, &ans, fp);
    //fun(10, alpha, 1, &ans);
    fprintf(fp,"%.*lf\n", NumOfSign, ans);

    return 0;
}
