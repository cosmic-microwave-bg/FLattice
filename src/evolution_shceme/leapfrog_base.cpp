#include "leapfrog_base.hpp"


double LeapFrogBase::laplacian ( const double* f, int i, int j, int k ) const
{	
    int ip1 = (i == N-1)?     0: i+1;
    int ip2 = (i >= N-2)? i-N+2: i+2;
    int im1 = (i ==   0)?   N-1: i-1;
    int im2 = (i <    2)? i+N-2: i-2;

    int jp1 = (j == N-1)?     0: j+1;
    int jp2 = (j >= N-2)? j-N+2: j+2;
    int jm1 = (j ==   0)?   N-1: j-1;
    int jm2 = (j <    2)? j+N-2: j-2;
    
    #if   DIMENSION == 1
        int idx = i;
        return (- f[ip2] + 16*f[ip1] - 30*f[idx] + 16*f[im1] - f[im2]) / (12*dx*dx);
    #elif DIMENSION == 2
        int idx = i*N + j;
        return ( (- f[ip2*N+j] + 16*f[ip1*N+j] - 30*f[idx] + 16*f[im1*N+j] - f[im2*N+j]) 
               + (- f[i*N+jp2] + 16*f[i*N+jp1] - 30*f[idx] + 16*f[i*N+jm1] - f[i*N+jm2]) ) / (12*dx*dx);
    #elif DIMENSION == 3
        int kp1 = (k == N-1)?     0: k+1;
        int kp2 = (k >= N-2)? k-N+2: k+2;
        int km1 = (k ==   0)?   N-1: k-1;
        int km2 = (k <    2)? k+N-2: k-2;

        int idx = (i*N + j)*N + k;
        return ( (- f[(ip2*N+j)*N+k] + 16*f[(ip1*N+j)*N+k] - 30*f[idx] + 16*f[(im1*N+j)*N+k] - f[(im2*N+j)*N+k])
               + (- f[(i*N+jp2)*N+k] + 16*f[(i*N+jp1)*N+k] - 30*f[idx] + 16*f[(i*N+jm1)*N+k] - f[(i*N+jm2)*N+k])
               + (- f[(i*N+j)*N+kp2] + 16*f[(i*N+j)*N+kp1] - 30*f[idx] + 16*f[(i*N+j)*N+km1] - f[(i*N+j)*N+km2]) ) / (12*dx*dx);
    #endif
}


void  LeapFrogBase::evolFields (double** f, double** df, const double h, double a)
{
    for( int n = 0; n < num_fields; ++n ){
        #if   DIMENSION == 1
        #pragma omp parallel for simd schedule(static) num_threads(num_threads)
        #elif DIMENSION >= 2
        #pragma omp parallel for schedule( static ) num_threads( num_threads )
        #endif
        for( int i = 0; i < N; ++i ){
            #if   DIMENSION == 1
                int idx = i;
                f[n][idx] += df[n][idx] * pow(a, A-D) * h*dt;
            #elif DIMENSION == 2
                #pragma omp simd
                for( int j = 0; j < N; ++j ){
                    int idx = i*N+j;
                    f[n][idx] += df[n][idx] * pow(a, A-D) * h*dt;
                }
            #elif DIMENSION == 3
                for( int j = 0; j < N; ++j ){
                    #pragma omp simd
                    for( int k = 0; k < N; ++k ){
                        int idx = (i*N+j)*N+k;
                        f[n][idx] += df[n][idx] * pow(a, A-D) * h*dt;
                    }
                }
            #endif
        }
    }
}


void LeapFrogBase::evolFieldDerivs(double** f, double** df, const double h, double a)
{
    for( int n = 0; n < num_fields; ++n ){
        #if DIMENSION == 1
        #pragma omp parallel for simd schedule(static) num_threads(num_threads)
        #elif DIMENSION >= 2
        #pragma omp parallel for schedule(static) num_threads(num_threads) 
        #endif
        for( int i = 0; i < N; ++i ){
            #if DIMENSION == 1
                int idx = i;
                df[n][idx] += ( laplacian(f[n], i)/pow(a, 2)  - dV(f, n, idx) ) * pow(a, A+D) * h*dt;
            #elif DIMENSION == 2
                #pragma omp simd
                for( int j = 0; j < N; ++j ){
                    int idx = i*N+j;
                    df[n][idx] += ( laplacian(f[n], i, j)/pow(a, 2)  - dV(f, n, idx) ) * pow(a, A+D) * h*dt;
                }
            #elif DIMENSION == 3
                for( int j = 0; j < N; ++j ){
                    #pragma omp simd
                    for( int k = 0; k < N; ++k ){
                        int idx = (i*N+j)*N+k;
                        df[n][idx] += ( laplacian(f[n], i, j, k)/pow(a, 2)  - dV(f, n, idx) ) * pow(a, A+D) * h*dt;
                    }
                }
            #endif
        }
    }
}