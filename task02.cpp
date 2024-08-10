#include "HW2.h"

using namespace std;

int main(int argc, char *argv[])
{
    double t1 = 0;
    double* Matrix, * b, * x, * y, * z,  * result, norm, norm_difference, n4;
    int n = 0, p = 0, r, s = 0, *d;
    char* filename = nullptr;
    if (!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n) == 1 && (n > 0)))
    {
        std::cout << "The entered numbers is not read" << std::endl;
        return -1;
    }
    p = atoi(argv[2]);
    r = atoi(argv[3]);
    s = atoi(argv[4]);
    if (argc == 6) filename = argv[5];

    Matrix = (double*)malloc(n * n * sizeof(double));
    if (!Matrix) { printf("Memory could not be allocated\n"); return 4; }
    b = (double*)malloc(n * sizeof(double));
    if (!b) { printf("Memory could not be allocated\n"); return 5; }
    d = (int*)malloc(n * sizeof(int));
    if (!b) { printf("Memory could not be allocated\n"); return 6; }
    x = (double*)malloc(n * sizeof(double));
    if (!x) { printf("Memory could not be allocated\n"); return 7; }
    y = (double*)malloc(n * sizeof(double));
    if (!y) { printf("Memory could not be allocated\n"); return 8; }
    z = (double*)malloc(n * sizeof(double));
    if (!z) { printf("Memory could not be allocated\n"); return 9; }
    result = (double*)malloc(n * sizeof(double));
    if (!result) { printf("Memory could not be allocated\n"); return 10; }
   
    //p > n ? p = n : p = p;
    if (Matrix_input(Matrix, b, n, s, filename)) {
        free(Matrix); free(d); free(b); free(x); free(y); free(z); free(result);
        std::cout << "Matrix is not input\n" << std::endl; 
        return -3;
    }
    printM(Matrix, r, n);
    printf("\n\n");
    n4=Norm(Matrix, n);
    if (method_Cholesky(n, Matrix, b, x, y, z,  d, p, t1, n4)!= 0) {
        free(Matrix); free(d); free(b); free(x); free(y); free(z); free(result);
        std::cout << "not corrected" << std::endl; return -4;
    }
    if (Matrix_input(Matrix, b, n, s, filename)) {
        free(Matrix); free(d); free(b); free(x); free(y); free(z); free(result);
        std::cout << "Matrix is not input\n" << std::endl; return -3;
    }
    Neviazka(Matrix, x, b, result, n, &norm, &norm_difference);
    printf("%s : residual = %e elapsed = %.2f n = %d p = %d r = %d s = %d\n",argv[0], norm, t1, n, p , r, s);
    free(Matrix); free(d); free(b); free(x); free(y); free(z); free(result);
}
