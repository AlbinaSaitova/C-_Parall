#include "HW2.h"
using namespace std;


int fscanMatrix(double* Matrix, int n, const char* filename) {
    int i, j;
    FILE* fin;
    if (!(fin = fopen(filename, "r"))) return -1;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (fscanf(fin, "%lg", &Matrix[i * n + j]) != 1) { return -2; }
        }
    }
    fclose(fin);
    return 0;
}
double function(int k, int n, int i, int j) {
    while (true) {
        switch (k) {
        case 1: {
            return double(n - (std::max(i, j)) + 1);
        }
        case 2: {
            return double(std::max(i, j));
        }
        case 3: {
            return double(abs(i - j));
        }
        case 4: {
            return double(1.0 / (i + j - 1));
        }
        }
    }
    return -4;
}
void vector_b(double* Matrix, double* b, int n) {
    int i, k;
    for (i = 0; i < n; i++) {
        b[i] = 0;
        for (k = 0; k < (n / 2) + (n % 2); k++) {

            b[i] += Matrix[i * n + 2 * k];

        }
    }
}
int Matrix_input(double* Matrix, double* b, int n, int k, const char* filename) {
    int i, j;
    if (k == 0) {
        if (fscanMatrix(Matrix, n, filename) != 0) {
            std::cout << "Error reading file" << std::endl; return -1;
        }
    }
    else if (k >= 1 && k <= 4) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                Matrix[i * n + j] = function(k, n, i + 1, j + 1);
            }
        }
    }
    else { std::cout << "The data is entered incorrectly" << std::endl; return 1; }
    vector_b(Matrix, b, n);
    return 0;
}

void printM(double* Matrix, int m, int n) {
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            printf("%10.3e ", Matrix[i * n + j]);
        }
        printf("\n");
    }
}
double Norm(double* a, int n){
    double n4=0.0;
    int j, i;
     for (j = 0; j < n; j++) {
		double s = 0;
		for (i = 0; i < n; i++) {
			s += a[i * n + j] * a[i * n + j];
		}
		n4 = ((s < n4) || n4 < 1.e-14 ? s : n4);
	}
	if(n4>1) n4=1;
	return n4;
    
}

double get_full_time()
{
    struct timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec / 1.e6;
}
int Sign(double a, double n4) {
    if (fabs(a) < n4*1.e-14)  return 0;
    if (a > 0)  return 1;
    else return -1;
}
void* t_f(void* args)
{
    Args* par = (Args*)args;
    int n = par->n;
    //int p = par->p;
    double n4 = par->n4;
    double* a = par->a;
    double* b = par->b;
    double* x = par->x;
    double* y = par->y;
    double* z = par->z;
    int* d = par->d;
    double t;
    int k, l;
    ///
    reduce_sum(par->p, 0, 0);
    for (l = 0; l < n; l++) {
         if (fabs(a[l * n + l]) <= n4*1.e-14) {
                par->err = 1;
                return nullptr;
            }
        if (par->i == 0) {
            t = a[l * n + l];
            for (k = 0; k <= l - 1; k++) {
                t -= a[k * n + l] * a[k * n + l] * d[k];
            }
            d[l] = Sign(t, n4);
            a[l* n + l] = sqrt(fabs(t));
        }
        reduce_sum(par->p, 0, 0);
       for (int j = par->i; j < n; j += par->p)
        {
           if (j < l + 1) continue;
            t = a[l * n + j];
           for (k = 0; k <= l - 1; k++) {
               t -= a[k * n + l] * d[k] * a[k * n + j];
            }
            if (fabs(a[l * n + l]) < n4*1.e-14 || fabs(d[l]) < n4*1.e-14) {
                 par->err = 2;
                return nullptr;
            }
            a[l * n + j] = t / (a[l * n + l] * d[l]);
        }
       reduce_sum(par->p, 0, 0);
    }
    reduce_sum(par->p, 0, 0);
    if (par->i == 0) {
        if (System_Solution1(a, z, b, n, n4) != 0){ par->err = 3;
            return nullptr;}
        if (System_Solution2(d, y, z, n, n4) != 0)  {par->err = 4;
        return nullptr;}
        if (System_Solution3(a, x, y, n, n4) != 0) { par->err = 5;
            return nullptr;}
    }

    return nullptr;
}

int method_Cholesky(int n, double* a, double* b, double* x, double* y, double* z, int* d, int p, double &t1, double n4){
        Args *par = new Args[p];
        pthread_t *tid = new pthread_t[p - 1];
        //pthread_t *tid = new pthread_t[p];
        for (int i = 0; i < p; i++)
        {
            par[i].a = a;
            par[i].b = b;
            par[i].d = d;
            par[i].x = x;
            par[i].y = y;
            par[i].z = z;
            par[i].i = i;
            par[i].n = n;
            par[i].p = p;
            par[i].n4=n4;
        }
        t1 = get_full_time();
         /*for (int i = 0; i < p; i++)
        {
           if( pthread_create(tid + i, 0, t_f, par + i)!=0){
            fprintf(stderr, "Cannot create thread %d\n", i);
            return 5;
               
            }
        }
        for (int i = 0; i < p ; i++)
        {
            if(pthread_join(tid[i], 0)!=0){
            fprintf(stderr, "Cannot wait thread %d\n", i);
            return 5;
            }
        }*/
        for (int i = 0; i < p - 1; i++)
        {
           if( pthread_create(tid + i, 0, t_f, par + i + 1)!=0){
            fprintf(stderr, "Cannot create thread %d\n", i);
             delete[] par;
             delete[] tid;
            return 5;
               
            }
        }
        t_f(par);
        for (int i = 0; i < p - 1; i++)
        {
            if(pthread_join(tid[i], 0)!=0){
            fprintf(stderr, "Cannot wait thread %d\n", i);
             delete[] par;
             delete[] tid;
            return 5;
            }
        }
        t1 = get_full_time() - t1;
        int res = par->err;
        delete[] par;
        delete[] tid;
        return res;
}
int System_Solution1(double* Matrix, double* x, double* b, int n, double n4) {
    int i, k;
    for (i = 0; i < n; i++) {
        x[i] = b[i];
        for (k = 0; k <= i - 1; k++) {
            x[i] -= Matrix[k * n + i] * x[k];
        }
        if (fabs(Matrix[i * n + i]) < n4*1.e-14) return 1;
        x[i] /= Matrix[i * n + i];
    }
    return 0;
}
int System_Solution2(int* d, double* y, double* z, int n, double n4) {
    int i;
    for (i = 0; i < n; i++) {
        if (fabs(d[i]) < n4*1.e-14) return 1;
        y[i] = z[i] / d[i];
    }
    return 0;
}
int System_Solution3(double* Matrix, double* x, double* y, int n, double n4) {
    int i, k;
    for (i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (k = i + 1; k < n; k++) {
            x[i] -= Matrix[i * n + k] * x[k];
        }
        if (fabs(Matrix[i * n + i]) < n4*1.e-14) return 1;
        x[i] /= Matrix[i * n + i];
    }
    return 0;
}

void Mult_M_V(double* a, int n, double* vector, double* result) {
    int i, j;
    for (i = 0; i < n; i++) {
        result[i] = 0.0;
        for (j = 0; j < n; j++) {
            result[i] += a[i * n + j] * vector[j];
        }
    }
}
void Neviazka(double* a, double* x, double* b, double* result, int  n, double* norm, double* norm_difference) {
    int i;
    double n1 = 0.0, n2 = 0.0, n3 = 0.0;
    Mult_M_V(a, n, x, result);
    for (i = 0; i < n; i++) {
        n1 += (result[i] - b[i]) * (result[i] - b[i]);
        n2 += b[i] * b[i];
    }
    *norm = sqrt(n1) / sqrt(n2);
    for (i = 0; i < n; i++) {
        n3 += (x[i] - ((i + 1) % 2)) * (x[i] - ((i + 1) % 2));
    }
    *norm_difference = sqrt(n3);
}
void reduce_sum(int total_threads, int *a, int k)
{

    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
    static int threads_in = 0;
    static int threads_out = 0;
    static int *pres = 0;
    int i;
    pthread_mutex_lock(&mutex);
    if (!pres)
    {
        pres = a;
    }
    else
    {
        for (i = 0; i < k; i++)
        {
            pres[i] += a[i];
        }
    }
    threads_in++;
    if (threads_in >= total_threads)
    {
        threads_out = 0;
        pthread_cond_broadcast(&condvar_in);
    }
    else
    {
        while (threads_in < total_threads)
        {
            pthread_cond_wait(&condvar_in, &mutex);
        }
    }
    if (pres != a)
    {
        for (i = 0; i < k; i++)
        {
            a[i] = pres[i];
        }
    }
    threads_out++;
    if (threads_out >= total_threads)
    {
        pres = 0;
        threads_in = 0;
        pthread_cond_broadcast(&condvar_out);
    }
    else
    {
        while (threads_out < total_threads)
        {
            pthread_cond_wait(&condvar_out, &mutex);
        }
    }
    pthread_mutex_unlock(&mutex);
}


