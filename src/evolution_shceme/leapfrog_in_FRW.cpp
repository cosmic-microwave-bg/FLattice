#include "leapfrog_in_FRW.hpp"



void LeapFrogRadMat::evolution(double** f, double** df, int output_step)
{
    if(_precision == 2)
    {
        evolFields(f, df, 0.5, a);
        t += 0.5*dt;
        a = pow(_coeff*t, _index);
        for(int i = 0; i < output_step-1; ++i){
            evolFieldDerivs(f, df, 1., a);
            t += 0.5*dt;
            a = pow(_coeff*t, _index);
            evolFields(f, df, 1., a);
            t += 0.5*dt;
            a = pow(_coeff*t, _index);
        }
        evolFields(f, df, 0.5, a);
        t += 0.5*dt;
        a = pow(_coeff*t, _index);
    }
    if(_precision == 4)
    {
        double t_f  = t;
        double t_df = t;
        for(int i = 0; i < output_step; ++i){
            evolFields(f, df, _C[0], a);
            for(int p = 0; p < 3; ++p){
                t_f +=  _C[p]*dt;
                a = pow(_coeff*t_f, _index);
                evolFieldDerivs(f, df, _D[p], a);
                t_df +=  _D[p]*dt;
                a = pow(_coeff*t_df, _index);
                evolFields(f, df, _C[p+1], a);
            }
            t_f += _C[3]*dt;
            a = pow(_coeff*t_f, _index);
        }
    }
}


void LeapFrogSelf::evolScale(double** f, const double h)
{
    a += da *h*dt;

    double C = 0;
    for( int n = 0; n < num_fields; ++n ){
        #if   DIMENSION == 1
        #pragma omp parallel for simd reduction(+:C) schedule(static) num_threads(num_threads)
        #elif DIMENSION >= 2
        #pragma omp parallel for reduction(+:C) schedule(static) num_threads(num_threads)
        #endif
        for( int i = 0; i < N; ++i ){
            #if   DIMENSION == 1
                int idx = i;
                C += 2*gradientEnergy(f[n], i)/(a*a) + 3*V(f, n, idx);
            #elif DIMENSION == 2
                #pragma omp simd reduction(+:C)
                for( int j = 0; j < N; ++j ){
                    int idx = i*N+j;
                    C += 2*gradientEnergy(f[n], i, j)/(a*a) + 3*V(f, n, idx);
                }
            #elif DIMENSION == 3
                for( int j = 0; j < N; ++j ){
                    #pragma omp simd reduction(+:C)
                    for( int k = 0; k < N; ++k ){
                        int idx = (i*N+j)*N+k;
                        C += 2*gradientEnergy(f[n], i, j, k)/(a*a) + 3*V(f, n, idx);
                    }
                }
            #endif
        }
    }
    C *= R*R*pow(a, 2*A+1)/3;
    C /= pow(N, DIMENSION);
    da = a/((A-2)*h*dt) * ( 1 - sqrt(1 - (A-2)*h*dt*(2*da + C*h*dt)/a) );
}

void LeapFrogSelf::calcDa(double** f, double** df)
{
    double rho = 0;
    for( int n = 0; n < num_fields; ++n ){
        int i = 0, j = 0, k = 0;
        #pragma omp parallel for reduction(+:rho) schedule(static) num_threads (num_threads)
        for( i = 0; i < N; ++i ){
            int idx = i;
        #if DIMENSION >= 2
            for( j = 0; j < N; ++j ){
                idx = i*N+j;
        #endif
        #if DIMENSION == 3
                for( k = 0; k < N; ++k ){
                    idx = (i*N+j)*N+k;
        #endif
                    rho += pow(df[n][idx]/pow(a,D), 2)/2 + gradientEnergy(f[n], i, j, k)/pow(a,2) + V(f, n, idx);
        #if DIMENSION == 3
                }
        #endif
        #if DIMENSION >= 2
            }
        #endif
        }
    }
    rho /= pow(N, DIMENSION);
    da = pow(a, A+1)*R*sqrt(rho/3);
}

void LeapFrogSelf::evolution( double** f, double** df, int output_step )
{
    calcDa(f, df);
    
    switch( _precision )
    {
        case 2:
            evolFields(f, df, 0.5, a);
            evolScale(f, 0.5);
            for( int i = 0; i < output_step-1; ++i )
            {
                evolFieldDerivs(f, df, 1.0, a);
                a += da *0.5*dt;
                evolFields(f, df, 1.0, a);
                evolScale(f, 0.5);
            }
            evolFields(f, df, 0.5, a);
            evolScale(f, 0.5);
            break;

        case 4:
            double a_df = a;
            for(int i = 0; i < output_step; ++i)
            {
                evolFields(f, df, _C[0], a);
                evolScale(f, _C[0]);
                for(int p = 0; p < 3; ++p){
                    evolFieldDerivs(f, df, _D[p], a_df);
                    a_df += da *_D[p]*dt;
                    evolFields(f, df, _C[p+1], a);
                    evolScale(f, _C[p+1]);
                }
            }
            break;
    }
}