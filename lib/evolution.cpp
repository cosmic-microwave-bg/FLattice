#include <cmath>
#include "parameter.hpp"
#include "utilities.hpp"
#include "evolution.hpp"

double LeapFrog::_a  = 1;
double LeapFrog::_da = 0;

LeapFrog::LeapFrog(): _dda()
{
    switch( precision )
    {
        case 2:
        case 4:
            break;
        default:
            Logout( "LeapFrog precision must be 2 or 4. \n" );
            exit(1);
    }
    
    switch( expansion )
    {
        case 0:
        case 1:
            break;
        case 2:
            _da  = 1;
            _dda = 0;
            break;
        case 3:
            _dda = 2;
            break;
        default:
            Logout( "Parameter 'expansion' must be 0 ~ 3 when you use class 'LeapFrog'. \n" );
            exit(1);
    }
}

void LeapFrog::evol_fields( double** f, double** df, double h )
{	
    for( int i = 0; i < num_fields; ++i ){
        #pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
            #if dim == 1
                int idx = j;
                f[i][idx] += df[i][idx] * h*dt;
            #elif dim == 2
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    f[i][idx] += df[i][idx] * h*dt;
                }
            #elif dim == 3
                for( int k = 0; k < N; ++k ){
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        f[i][idx] += df[i][idx] * h*dt;
                    }
                }
            #endif
        }
    }
}

void LeapFrog::evol_field_derivs( double** f, double** df, Field* field, double h )
{
    for( int i = 0; i < num_fields; ++i ){
        #ifdef SPHERICAL_SYM
            df[i][N-1] -= ( field->laplacian(f[i], df[i], N-1, 0, 0) + f[i][N-1]/2) * h*dt;
            #pragma omp parallel for schedule( static ) num_threads(num_threads)
            for( int j = 0; j < N-1; ++j ){
		#else
            #pragma omp parallel for schedule( static ) num_threads( num_threads )
            for( int j = 0; j < N; ++j ){
        #endif
            #if dim == 1
                int idx = j;
                df[i][idx] += ( field->laplacian(f[i], j) - field->dV(f, i, idx) ) * h*dt;
            #elif dim == 2
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    df[i][idx] += ( field->laplacian(f[i], j, k) - field->dV(f, i, idx) ) * h*dt;
                }
            #elif dim == 3
                for( int k = 0; k < N; ++k ){
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        df[i][idx] += ( field->laplacian(f[i], j, k, l) - field->dV(f, i, idx) ) * h*dt;
                    }
                }
            #endif
        }
    }
}

void LeapFrog::evol_field_derivs_expansion( double** f, double** df, Field* field, double h )
{
    for( int i = 0; i < num_fields; ++i ){
        #pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
            #if dim == 1
                int idx = j;
                df[i][idx] += ( field->laplacian(f[i], j, 0, 0) + _dda*f[i][idx]/_a - pow(_a,3)*field->adV(f, i, idx, _a) ) * h*dt;
            #elif dim == 2
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    //df[i][idx] += ( field->laplacian(f[i], j, k, 0) + _dda*f[i][idx]/_a - pow(_a,3)* (field->*(field->adV[i]))(f, _a, i, idx) ) * h*dt;
                    df[i][idx] += ( field->laplacian(f[i], j, k, 0) + _dda*f[i][idx]/_a - pow(_a,3)*field->adV(f, i, idx, _a) ) * h*dt;
                }
            #elif dim == 3
                for( int k = 0; k < N; ++k ){
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        df[i][idx] += ( field->laplacian(f[i], j, k, l) + _dda*f[i][idx]/_a - pow(_a,3)*field->adV(f, i, idx, _a) ) * h*dt;
                    }
                }
            #endif
        }
    }
}


void LeapFrog::evol_scale( Field* field, double** f, double h )
{
    _a += _da *h*dt;
    
    double C = 0;
    C += field->gradient_energy(f[0])/(3*_a) + pow(_a,3)*field->potential_energy( f[0], _a );
    for( int i = 0; i < dim; ++i ) C /= N;

    _dda = - 2./(h*dt) * ( _da + _a/(h*dt) * ( 1 - sqrt( 1 + ( 2*h*dt*_da + pow(h*dt,2)*C )/_a ) ) );
}


void LeapFrog::evolution( Field* field, double** f, double** df )
{
    if( expansion != 0 )
    {
        Logout( "Error: Parameter expansion must be 0 when you use member-function 'evolution'. \n" );
        Logout( "Use expansion = 1 ~ 3 and member-function 'evolution_expansion' when you simulate the expanding universe. \n" );
        exit(1);
    }
    
    switch( precision )
    {
        case 2:
            evol_fields( f, df, 0.5 );
            for( int i = 0; i < output_step; ++i )
            {
                evol_field_derivs( f, df, field, 1.0 );
                if( i == output_step - 1 ) evol_fields( f, df, 0.5 );
                else evol_fields( f, df, 1.0 );
            }
            break;

        case 4:
            const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
            const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
            for( int i = 0; i < output_step; ++i )
            {
                evol_fields( f, df, C[0] );
                for( int p = 0; p < 3; ++p )
                {
                    evol_field_derivs( f, df, field, D[p] );
                    evol_fields( f, df, C[p+1] );
                }
            }
            break;
    }
}

void LeapFrog::evolution_expansion( Field* field, double** f, double** df, double t )
{   
    switch ( expansion )
    {   
        /* // The code is not well optimized when include member-function 'evol_field_derivs'
        case 0: // no expansion
            switch( precision )
            {
                case 2:
                    evol_fields( f, df, 0.5 );
                    for( int i = 0; i < output_step; ++i )
                    {
                        evol_field_derivs( f, df, field, 1.0 );
                        if( i == output_step - 1 ) evol_fields( f, df, 0.5 );
                        else evol_fields( f, df, 1.0 );
                    }
                    break;

                case 4:
                    const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
                    const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
                    for( int i = 0; i < output_step; ++i )
                    {
                        evol_fields( f, df, C[0] );
                        for( int p = 0; p < 3; ++p ){
                            evol_field_derivs( f, df, field, D[p] );
                            evol_fields( f, df, C[p+1] );
                        }
                    }
                    break;
            }
            break;
        */
        case 1: // self-consistent evolution
            _da = sqrt( ( field->gradient_energy(f[0]) + 2*pow(_a,4)*field->potential_energy(f[0], _a) )/6 );
            for( int i = 0; i < dim; ++i ) _da /= sqrt(N);

            switch( precision )
            {
                case 2:
                    evol_fields( f, df, 0.5 );
                    evol_scale( field, f, 0.5 );
                    for( int i = 0; i < output_step; ++i )
                    {
                        evol_field_derivs_expansion( f, df, field, 1.0 );
                        evol_scale_derivs( 1.0 );
                        if( i == output_step - 1 )
                        {
                            evol_fields( f, df, 0.5 );
                            evol_scale( field, f, 0.5 );
                        }
                        else
                        {
                            evol_fields( f, df, 1.0 );
                            evol_scale( field, f, 1.0 );
                        }
                    }
                    break;

                case 4:
                    const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
                    const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
                    for( int i = 0; i < output_step; ++i )
                    {
                        evol_fields( f, df, C[0] );
                        evol_scale( field, f, C[0] );
                        for( int p = 0; p < 3; ++p ){
                            evol_field_derivs_expansion( f, df, field, D[p] );
                            evol_scale_derivs( D[p] );
                            evol_fields( f, df, C[p+1] );
                            evol_scale( field, f, C[p+1] );
                        }
                    }
                    break;
            }
            break;

        case 2: // radiation dominated universe
            switch( precision )
            {
                case 2:
                    evol_fields( f, df, 0.5 );
                    _a +=  0.5*dt;
                    for( int i = 0; i < output_step; ++i ){
                        evol_field_derivs_expansion( f, df, field, 1. );
                        if( i == output_step-1 )
                        {
                            evol_fields( f, df, 0.5 );
                            _a +=  0.5*dt;
                        }
                        else
                        {
                            evol_fields( f, df, 1. );
                            _a +=  dt;
                        }
                    }
                    break;

                case 4:
                    const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
                    const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
                    for( int i = 0; i < output_step; ++i ){
                        evol_fields( f, df, C[0] );
                        for( int p = 0; p < 3; ++p ) {
                            _a +=  C[p]*dt;
                            evol_field_derivs_expansion( f, df, field, D[p] );
                            evol_fields( f, df, C[p+1] );
                        }
                        _a +=  C[3]*dt;
                    }
                    break;
            }
            break;

        case 3: // matter dominated universe
            switch( precision )
            {
                case 2:
                    evol_fields( f, df, 0.5 );
                    t +=  0.5*dt;
                    _a = pow( t, 2 );
                    for( int i = 0; i < output_step; ++i ){
                        evol_field_derivs_expansion( f, df, field, 1.0 );
                        if( i == output_step - 1 )
                        {
                            evol_fields( f, df, 0.5 );
                            t +=  0.5*dt;
                        }
                        else
                        {
                            evol_fields( f, df, 1.0 );
                            t +=  dt;
                        }
                        _a = pow( t, 2 );
                    }
                    _da = 2*t;
                    break;

                case 4:
                    const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
                    const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
                    for( int i = 0; i < output_step; ++i ){
                        evol_fields( f, df, C[0] );
                        for( int p = 0; p < 3; ++p ) {
                            t +=  C[p]*dt;
                            _a = pow( t, 2 );
                            evol_field_derivs_expansion( f, df, field, D[p] );
                            evol_fields( f, df, C[p+1] );
                        }
                        t +=  C[3]*dt;
                        _a = pow( t, 2 );
                    }
                    _da = 2*t;
                    break;
            }
            break;
        
        default:
            Logout( "Error: Parameter 'expansion' must be 0 ~ 3. \n" );
            exit(1);
    }
}