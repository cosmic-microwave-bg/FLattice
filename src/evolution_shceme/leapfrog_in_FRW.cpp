#include "leapfrog_in_FRW.hpp"



void LeapFrogRadMat::evolution( double** f, double** df, int output_step )
{
    // a ~ t in radiation dominated univ
    // a ~ t^2 in matter dominated univ
    int exp = (expansion==Expansion::rad_domination)? 1: 2;
    switch( _precision )
    {
        case 2:
            evolFields( f, df, 0.5 );
            t +=  0.5*dt;
            a =  pow(t, exp);
            for( int i = 0; i < output_step; ++i ){
                evolFieldDerivs( f, df, 1., a);
                if( i == output_step-1 ){
                    evolFields( f, df, 0.5 );
                    t +=  0.5*dt;
                }else{
                    evolFields( f, df, 1. );
                    t +=  dt;
                }
                a = pow(t,exp);
            }
            break;
        case 4:
            const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
            const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };

            for( int i = 0; i < output_step; ++i ){
                evolFields( f, df, C[0] );
                for( int p = 0; p < 3; ++p ) {
                    t +=  C[p]*dt;
                    a = pow(t, exp);
                    evolFieldDerivs( f, df, D[p], a );
                    evolFields( f, df, C[p+1] );
                }
                t +=  C[3]*dt;
                a = pow( t, exp );
            }
            break;
    }
    if( expansion == Expansion::mat_domination ) da = 2*t;
}


void LeapFrogSelf::evolScale( double** f, const double h )
{
    a += da *h*dt;

    double C = 0;
    for( int n = 0; n < num_fields; ++n ){
        #pragma omp parallel for reduction(+:C) schedule( static ) num_threads( num_threads )
        for( int i = 0; i < N; ++i ){
            #if   DIMENSION == 1
                int idx = i;
                C += gradientEnergy(f[n], i)/(3+a) + pow(a,3)*_model->V(f, n, idx, a);
            #elif DIMENSION == 2
                for( int j = 0; j < N; ++j ){
                    int idx = i*N+j;
                    C += gradientEnergy(f[n], i, j)/(3*a) + pow(a,3)*_model->V(f, n, idx, a);
                }
            #elif DIMENSION == 3
                for( int j = 0; j < N; ++j ){
                    for( int k = 0; k < N; ++k ){
                        int idx = (i*N+j)*N+k;
                        C += gradientEnergy(f[n], i, j, k)/(3*a) + pow(a,3)*_model->V(f, n, idx, a);
                    }
                }
            #endif
        }
    }
    for( int d = 0; d < DIMENSION; ++d ) C /= N;
    dda = - 2./(h*dt) * ( da + a/(h*dt) * (1 - sqrt( 1 + (2*h*dt*da + C*pow(h*dt,2))/a )) );
}


void LeapFrogSelf::evolScaleDerivs(const double h)
{
    da += dda * h*dt;
}


void LeapFrogSelf::evolution( double** f, double** df, int output_step )
{
    switch( _precision )
    {
        case 2:
            evolFields( f, df, 0.5 );
            evolScale( f, 0.5 );
            for( int i = 0; i < output_step; ++i )
            {
                evolFieldDerivs( f, df, 1.0, a );
                evolScaleDerivs( 1.0 );
                if( i == output_step - 1 ){
                    evolFields( f, df, 0.5 );
                    evolScale( f, 0.5 );
                }else{
                    evolFields( f, df, 1.0 );
                    evolScale( f, 1.0 );
                }
            }
            break;

        case 4:
            const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
            const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };

            for( int i = 0; i < output_step; ++i )
            {
                evolFields( f, df, C[0] );
                evolScale( f, C[0] );
                for( int p = 0; p < 3; ++p ){
                    evolFieldDerivs( f, df, D[p], a );
                    evolScaleDerivs( D[p] );
                    evolFields( f, df, C[p+1] );
                    evolScale( f, C[p+1] );
                }
            }
            break;
    }
}