//
// Created by osboxes on 12/2/21.
//


#ifndef UNTITLED1_INTEGR_H
#define UNTITLED1_INTEGR_H

typedef int(*function)(double x, double *ans);
typedef int(*integral)(function fx, double a, double b, double eps, unsigned long long steps, double *ans);


int countSteps(integral I, function fx, double a, double b, int m, double eps, unsigned long long *ans);
int Integral_midPoint(function fx, double a, double b, double eps, unsigned long long steps, double *ans);
int Integral_Gauss(function fx, double a, double b, double eps, unsigned long long steps, double *ans);
int Integral_Simpson(function fx, double a, double b, double eps, unsigned long long steps, double *ans);

#endif //UNTITLED1_INTEGR_H
