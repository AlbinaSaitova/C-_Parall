#include <stdio.h>
#include <iostream>
#include <math.h>
#include <limits>
#include <pthread.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/sysinfo.h>
#include <sys/resource.h>

struct Args
{
    int n = 0;
    int i = 0;
    int p = 0;
    char *filename;
    double *a;
    double *b;
    double* x;
    double* y;
    double* z;
    int *d;
    int err = 0;
    double t = 0;
    double n4 = 0;
    double r1 = 0;
    double r2 = 0;
    ~Args()
    {
        a = nullptr;
        b = nullptr;
        d = nullptr;
        x = nullptr;
        y = nullptr;
        z = nullptr;
    }
};
int Matrix_input(double* Matrix, double* b, int n, int k, const char* filename);
void reduce_sum(int total_threads, int *a, int k);
void printM(double* Matrix, int m, int n);
double Norm(double* a, int n);
int method_Cholesky(int n, double* a, double* b, double* x, double* y, double* z, int* d, int p, double &t1, double n4);
int System_Solution1(double* Matrix, double* x, double* b, int n, double n4);
int System_Solution2(int* d, double* y, double* z, int n, double n4);
int System_Solution3(double* Matrix, double* x, double* b, int n, double n4);
void Neviazka(double* a, double* x, double* b, double* result, int  n, double* norm, double* norm_difference);
double get_full_time();



/*double get_full_time()
{
    struct timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec / 1.e6;
}*/


