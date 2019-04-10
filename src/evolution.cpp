#include <cmath>
#include "evolution.hpp"
#include "utilities.hpp"



double Evolution::laplacian ( const double* f, int j, int k, int l ) const
{	
    int jp1 = (j == N-1)?     0: j+1;
	int jp2 = (j >= N-2)? j-N+2: j+2;
	#ifdef SPHERICAL_SYM
	if( DIMENSION != 1 ){
		Logout( "Error: Dimension must be 1 when you define SPHERICAL_SYM with ABC. \n" );
		exit(1);
	}
	int jm1 = abs( j-1 );
	int jm2 = abs( j-2 );
	#else
	int jm1 = (j == 0)?   N-1: j-1;
	int jm2 = (j <  2)? j+N-2: j-2;
	#endif
	
	#if DIMENSION == 1
	    int idx = j;
		#ifdef SPHERICAL_SYM
            if( idx == 0 ){ // r = 0では1/r*df/dr = 0 より別に処理
                return (- f[jp2] + 16*f[jp1] - 30*f[idx] + 16*f[jm1] - f[jm2]) / (12*dx*dx);
            }else if( idx == N-2 ){ //境界1つ手前は中央差分2次で計算
                return (f[jp1] - 2*f[idx] + f[jm1]) / (dx*dx) + 2*gradient(f, 0, i, idx, 0, 0) / (idx*dx);
            }else if( idx == N-1 ){ //境界では後退差分2次で計算、時間発展はABC
                return (df[jm2] - 4*df[jm1] + 3*df[idx]) / (2*dx) + df[idx] / (idx*dx);
            }else{
                return (- f[jp2] + 16*f[jp1] - 30*f[idx] + 16*f[jm1] - f[jm2]) / (12*dx*dx) + 2*gradient(f, 0, i, idx, 0, 0) / (idx*dx);
            }
		#else
		    return (- f[jp2] + 16*f[jp1] - 30*f[idx] + 16*f[jm1] - f[jm2]) / (12*dx*dx);
		#endif
	#elif DIMENSION == 2
		int kp1 = (k == N-1)?     0: k+1;
		int kp2 = (k >= N-2)? k-N+2: k+2;
		int km1 = (k ==   0)?   N-1: k-1;
		int km2 = (k <    2)? k+N-2: k-2;

        int idx = j*N + k;
        return ( (- f[jp2*N+k] + 16*f[jp1*N+k] - 30*f[idx] + 16*f[jm1*N+k] - f[jm2*N+k]) 
               + (- f[j*N+kp2] + 16*f[j*N+kp1] - 30*f[idx] + 16*f[j*N+km1] - f[j*N+km2]) ) / (12*dx*dx);
	#elif DIMENSION == 3
		int kp1 = (k == N-1)?     0: k+1;
		int kp2 = (k >= N-2)? k-N+2: k+2;
		int km1 = (k ==   0)?   N-1: k-1;
		int km2 = (k <    2)? k+N-2: k-2;
		
		int lp1 = (l == N-1)?     0: l+1;
		int lp2 = (l >= N-2)? l-N+2: l+2;
		int lm1 = (l ==   0)?   N-1: l-1;
		int lm2 = (l <    2)? l+N-2: l-2;

        int idx = (j*N + k)*N + l;
        return ( (- f[(jp2*N+k)*N+l] + 16*f[(jp1*N+k)*N+l] - 30*f[idx] + 16*f[(jm1*N+k)*N+l] - f[(jm2*N+k)*N+l])
               + (- f[(j*N+kp2)*N+l] + 16*f[(j*N+kp1)*N+l] - 30*f[idx] + 16*f[(j*N+km1)*N+l] - f[(j*N+km2)*N+l])
               + (- f[(j*N+k)*N+lp2] + 16*f[(j*N+k)*N+lp1] - 30*f[idx] + 16*f[(j*N+k)*N+lm1] - f[(j*N+k)*N+lm2]) ) / (12*dx*dx);
	#endif
}


void   Evolution::evol_fields ( double** f, double** df, double h )
{	
    for( int n = 0; n < num_fields; ++n ){
        #pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int i = 0; i < N; ++i ){
            #if   DIMENSION == 1
                int idx = i;
                f[n][idx] += df[n][idx] * h*dt;
            #elif DIMENSION == 2
                for( int j = 0; j < N; ++j ){
                    int idx = i*N+j;
                    f[n][idx] += df[n][idx] * h*dt;
                }
            #elif DIMENSION == 3
                for( int j = 0; j < N; ++j ){
                    for( int k = 0; k < N; ++k ){
                        int idx = (i*N+j)*N+k;
                        f[n][idx] += df[n][idx] * h*dt;
                    }
                }
            #endif
        }
    }
}

void   LeapFrog::evol_field_derivs( double** f, double** df, double h )
{
    for( int n = 0; n < num_fields; ++n ){
        #ifdef SPHERICAL_SYM
            df[n][N-1] -= ( laplacian(f[n], df[n], N-1, 0, 0) + f[n][N-1]/2) * h*dt;
            #pragma omp parallel for schedule( static ) num_threads(num_threads)
            for( int i = 0; i < N-1; ++i ){
		#else
            #pragma omp parallel for schedule( static ) num_threads( num_threads )
            for( int i = 0; i < N; ++i ){
        #endif
            #if DIMENSION == 1
                int idx = i;
                df[n][idx] += ( laplacian(f[n], i) - _model->dV(f, n, idx) ) * h*dt;
            #elif DIMENSION == 2
                for( int j = 0; j < N; ++j ){
                    int idx = i*N+j;
                    df[n][idx] += ( laplacian(f[n], i, j) - _model->dV(f, n, idx) ) * h*dt;
                }
            #elif DIMENSION == 3
                for( int j = 0; j < N; ++j ){
                    for( int k = 0; k < N; ++k ){
                        int idx = (i*N+j)*N+k;
                        df[n][idx] += ( laplacian(f[n], i, j, k) - _model->dV(f, n, idx) ) * h*dt;
                    }
                }
            #endif
        }
    }
}


void LeapFrog::evolution( double** f, double** df )
{
    if( expansion != Expansion::no_expansion )
    {
        Logout( "Error: Parameter expansion must be 'no_expansion' when you use member-function 'evolution'. \n" );
        exit(1);
    }
    
    switch( precision )
    {
        case 2:
            evol_fields( f, df, 0.5 );
            for( int i = 0; i < output_step; ++i )
            {
                evol_field_derivs( f, df, 1.0 );
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
                    evol_field_derivs( f, df, D[p] );
                    evol_fields( f, df, C[p+1] );
                }
            }
            break;
    }
}


double LeapFrogExpansion::_a   = 1;
double LeapFrogExpansion::_da  = 1;
double LeapFrogExpansion::_dda = 0;


LeapFrogExpansion::LeapFrogExpansion( std::shared_ptr<Model> model ): Evolution(model)
{
    switch( expansion )
    {
        case Expansion::self_consist:
            break;
        case Expansion::rad_domination:
            _dda = 0;
            break;
        case Expansion::mat_domination:
            _dda = 2;
            break;
        default:
            std::cerr << "Parameter 'expansion' must be enum class ... when you use class 'LeapFrog'." << std::endl;
            exit(1);
    }
}


void LeapFrogExpansion::evol_field_derivs( double** f, double** df, double h )
{
    for( int i = 0; i < num_fields; ++i ){
        #pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
            #if DIMENSION == 1
                int idx = j;
                df[i][idx] += ( laplacian(f[i], j) + _dda*f[i][idx]/_a - pow(_a,3)*_model->adV(f, i, idx, _a) ) * h*dt;
            #elif DIMENSION == 2
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    df[i][idx] += ( laplacian(f[i], j, k) + _dda*f[i][idx]/_a - pow(_a,3)*_model->adV(f, i, idx, _a) ) * h*dt;
                }
            #elif DIMENSION == 3
                for( int k = 0; k < N; ++k ){
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        df[i][idx] += ( laplacian(f[i], j, k, l) + _dda*f[i][idx]/_a - pow(_a,3)*_model->adV(f, i, idx, _a) ) * h*dt;
                    }
                }
            #endif
        }
    }
}


void LeapFrogExpansion::evol_scale( double** f, double h )
{
    _a += _da *h*dt;
    
    double C = 0;
    C += _model->gradient_energy(f[0])/(3*_a) + pow(_a,3)*_model->potential_energy( f[0], _a );
    for( int i = 0; i < DIMENSION; ++i ) C /= N;

    _dda = - 2./(h*dt) * ( _da + _a/(h*dt) * ( 1 - sqrt( 1 + ( 2*h*dt*_da + pow(h*dt,2)*C )/_a ) ) );
}


void LeapFrogExpansion::evolution( double** f, double** df, double t )
{   
    switch ( expansion )
    {   
        case Expansion::self_consist: // self-consistent evolution
            _da = sqrt( ( _model->gradient_energy(f[0]) + 2*pow(_a,4)*_model->potential_energy(f[0], _a) )/6 );
            for( int i = 0; i < DIMENSION; ++i ) _da /= sqrt(N);

            switch( precision )
            {
                case 2:
                    evol_fields( f, df, 0.5 );
                    evol_scale( f, 0.5 );
                    for( int i = 0; i < output_step; ++i )
                    {
                        evol_field_derivs( f, df, 1.0 );
                        evol_scale_derivs( 1.0 );
                        if( i == output_step - 1 )
                        {
                            evol_fields( f, df, 0.5 );
                            evol_scale( f, 0.5 );
                        }
                        else
                        {
                            evol_fields( f, df, 1.0 );
                            evol_scale( f, 1.0 );
                        }
                    }
                    break;

                case 4:
                    const double C[4] = { +0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844 };
                    const double D[3] = { +1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688 };
                    for( int i = 0; i < output_step; ++i )
                    {
                        evol_fields( f, df, C[0] );
                        evol_scale( f, C[0] );
                        for( int p = 0; p < 3; ++p ){
                            evol_field_derivs( f, df, D[p] );
                            evol_scale_derivs( D[p] );
                            evol_fields( f, df, C[p+1] );
                            evol_scale( f, C[p+1] );
                        }
                    }
                    break;
            }
            break;

        case Expansion::rad_domination: // radiation dominated universe
            switch( precision )
            {
                case 2:
                    evol_fields( f, df, 0.5 );
                    _a +=  0.5*dt;
                    for( int i = 0; i < output_step; ++i ){
                        evol_field_derivs( f, df, 1. );
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
                            evol_field_derivs( f, df, D[p] );
                            evol_fields( f, df, C[p+1] );
                        }
                        _a +=  C[3]*dt;
                    }
                    break;
            }
            break;

        case Expansion::mat_domination: // matter dominated universe
            switch( precision )
            {
                case 2:
                    evol_fields( f, df, 0.5 );
                    t +=  0.5*dt;
                    _a = pow( t, 2 );
                    for( int i = 0; i < output_step; ++i ){
                        evol_field_derivs( f, df, 1.0 );
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
                            evol_field_derivs( f, df, D[p] );
                            evol_fields( f, df, C[p+1] );
                        }
                        t +=  C[3]*dt;
                        _a = pow( t, 2 );
                    }
                    _da = 2*t;
                    break;
            }
            break;
    }
}