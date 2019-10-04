/**
 * @file  evolution.cpp
 * @brief See the coresponding header file for details.
 */
#include "leapfrog_in_Minkowsky.hpp"



void LeapFrogInMinkowsky::evolution( double** f, double** df, int output_step )
{
    switch( _precision )
    {
        case 2:
            evolFields( f, df, 0.5 );
            for( int i = 0; i < output_step; ++i ){
                evolFieldDerivs( f, df, 1.0 );
                if( i == output_step-1 ) evolFields( f, df, 0.5 );
                else evolFields( f, df, 1.0 );
            }
            break;

        case 4:
            for( int i = 0; i < output_step; ++i ){
                evolFields(f, df, _C[0]);
                for(int p = 0; p < 3; ++p){
                    evolFieldDerivs(f, df, _D[p]);
                    evolFields(f, df, _C[p+1]);
                }
            }
            break;
    }
}


double LeapFrogWithABC::sphericalSymLaplacian( double* f, int i )
{
    int ip1 = (i == N-1)?     0: i+1;
    int ip2 = (i >= N-2)? i-N+2: i+2;
    int im1 = abs( i-1 );  // for center 
	int im2 = abs( i-2 );  // for center
    
    int idx = i;
    double r = i*dx;

    if( idx == 0 )  // 原点では gradient term は 0
        return (- f[ip2] + 16*f[ip1] - 30*f[idx] + 16*f[im1] - f[im2] )/(12*dx*dx);
    if( idx == 1 )
        return (- f[ip2] + 16*f[ip1] - 30*f[idx] + 16*f[im1] - f[im2] )/(12*dx*dx) + 2/r*(f[ip1] - f[im1])/(2*dx);
    if( idx == 1 or idx == N-2 )  // 原点隣 or 境界1個手前は中央差分2次
        return (f[ip1] - 2*f[idx] + f[im1])/(dx*dx) + 2/r*(f[ip1] - f[im1])/(2*dx);
    if( idx == N-1 )  //境界では後退差分2次 (ABC用にgradientのみ)
        return (f[im2] - 4*f[im1] + 3*f[idx]) / (2*dx);

    return (- f[ip2] + 16*f[ip1] - 30*f[idx] + 16*f[im1] - f[im2])/(12*dx*dx)
     + 2/r*(- f[ip2] +  8*f[ip1]             -  8*f[im1] + f[im2])/(12*dx);
}

void LeapFrogWithABC::evolFieldDerivsWithABC( double **f, double **df, const double h )
{
    for( int n = 0; n < num_fields; ++n ){
        df[n][N-1] = - ( (df[n][N-3]-4*df[n][N-2] + 2*meff2*f[n][N-1]*dx)*h*dt + df[n][N-1]*(3*h*dt-4*dx) );
        #pragma omp parallel for simd schedule(static) num_threads(num_threads)
        for( int i = 0; i < N-1; ++i ){
            int idx = i;
            df[n][idx] += ( sphericalSymLaplacian(f[n], i) - dV(f, n, idx) ) * h*dt;
        }
        df[n][N-1] = ( df[n][N-1] - (df[n][N-3] - 4*df[n][N-2])*h*dt ) / (4*dx+3*h*dt);
        //df[n][N-1] -= ( sphericalSymLaplacian(df[n], N-1) + meff2*f[n][N-1]/2 + df[n][N-1]/L ) * h*dt;
    }
}

void LeapFrogWithABC::evolution( double** f, double** df, int loop )
{
    switch( _precision )
    {
        case 2:
            evolFields(f, df, 0.5);
            for(int i = 0; i < output_step; ++i){
                evolFieldDerivsWithABC(f, df, 1. );
                if( i == output_step-1 ) evolFields(f, df, 0.5);
                else evolFields(f, df, 1.0);
            }
            break;

        case 4:
            for(int i = 0; i < output_step; ++i){
                evolFields(f, df, _C[0]);
                for(int p = 0; p < 3; ++p){
                    evolFieldDerivsWithABC(f, df, _D[p]);
                    evolFields(f, df, _C[p+1]);
                }
            }
            break;
    }
}