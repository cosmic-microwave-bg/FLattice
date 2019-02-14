#include <cmath>

#include "calculation.hpp"
#include "parameter.hpp"


Calculation::Calculation( int num_fields ): _num_fields(num_fields)
{
    _total_average = 0;
    _average = new double [_num_fields]();
    _variance = new double [_num_fields]();

    value = new double* [_num_fields];
	switch( dim ){
		case 1:
			value[0] = new double [num_fields*N]();
            total_value = new double [N]();
			for( int i = 0; i < num_fields; ++i ) value[i] = value[0] + i*N;
			break;
		case 2:
			value[0] = new double [num_fields*N*N]();
            total_value = new double [N*N]();
			for( int i = 0; i < num_fields; ++i ) value[i] = value[0] + i*N*N;
			break;
		case 3:
			value[0] = new double [num_fields*N*N*N]();
            total_value = new double [N*N*N]();
			for( int i = 0; i < num_fields; ++i ) value[i] = value[0] + i*N*N*N;
			break;
	}
}


Energy::Energy( Field* field, LeapFrog* leapfrog, double** f, double** df, int num_fields, double dx ): _dx(dx), Calculation(num_fields)
{
    double a = leapfrog->a();
    double da = leapfrog->da();

    if( expansion )
    {
        for( int i = 0; i < _num_fields; ++i ){
            #pragma omp parallel for schedule( static ) num_threads ( num_threads )
            for( int j = 0; j < N; ++j ){
                switch( dim ){
                    case 1:
                        int idx = j;
                        value[i][idx] = pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,4)) + field->aV(f, a, i ,dx);
                        _average[i] += value[i][idx];
                        break;
                    case 2:
                        for( int k = 0; k < N; ++k ){
                            int idx = j*N + k;
                            value[i][idx] = pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,4)) + field->aV(f, a, i ,dx);
                            _average[i] += value[i][idx];
                        }
                        break;
                    case 3:
                        for( int k = 0; k < N; ++k ){
                            for( int l = 0; l < N; ++l ){
                                int idx = (j*N + k)*N + l;
                                value[i][idx] = pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,4)) + field->aV(f, a, i ,dx);
                                _average[i] += value[i][idx];
                            }
                        }
                        break;
                }
            }
        }
        for( int i = 0; i < _num_fields; ++i ){
            _average[i] += gradient_energy(f[i]) / (2*pow(a,4));
            for( int j = 0; j < dim; ++j ) _average[i] /= N;
            _total_average += _average[i];
        }
    }
    else
    {
        for( int i = 0; i < _num_fields; ++i ){
            #pragma omp parallel for schedule( static ) num_threads ( num_threads )
            for( int j = 0; j < N; ++j ){
                switch( dim ){
                    case 1:
                        int idx = j;
                        value[i][idx] = pow(df[i][idx],2) + field->V(f, i ,dx);
                        _average[i] += value[i][idx];
                        break;
                    case 2:
                        for( int k = 0; k < N; ++k ){
                            int idx = j*N + k;
                            value[i][idx] = pow(df[i][idx],2) + field->V(f, i ,dx);
                            _average[i] += value[i][idx];
                        }
                        break;
                    case 3:
                        for( int k = 0; k < N; ++k ){
                            for( int l = 0; l < N; ++l ){
                                int idx = ( j*N + k)*N + l;
                                value[i][idx] = pow(df[i][idx],2) + field->V(f, i ,dx);
                                _average[i] += value[i][idx];
                            }
                        }
                        break;
                }
            }
        }
        for( int i = 0; i < _num_fields; ++i ){
            _average[i] += gradient_energy( f[i] )/2;
            for( int j = 0; j < dim; ++j ) _average[i] /= N;
            _total_average += _average[i];
        }
    }
}


double Energy::gradient_energy( double* f )
{
    double gradient_energy = 0;
    
    #pragma omp parallel for schedule( static ) num_threads ( num_threads )
    for( int j = 0; j < N; ++j ){
        int jp1 = ( j ==N-1 ) ? 0 : j+1;
        int jp2 = ( j>=N-2 ) ? j-N+2 : j+2;
        #ifdef SPHERICAL_SYM
        int jm1 = abs( j-1 );
        int jm2 = abs( j-2 );
        #else
        int jm1 = ( j == 0 ) ? N-1 : j-1;
        int jm2 = ( j < 2 ) ? j+N-2 : j-2;
        #endif
        #if dim == 1
		    #ifdef SPHERICAL_SYM
            int idx = j;
            if( idx == 0 ){
                grad[0] = 0;
            }else if( idx == N-2 ){
                grad[0] = ( f[jp1] - f[jm1] ) / (2*_dx);
            }else if( idx == N-1 ){
                grad[0] = ( f[jm2] - 4*f[jm2] + 3*f[idx]  ) / (2*_dx);
            }else{
                grad[0] = ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*_dx);
            }
            #else
            gradient_energy +=  pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*_dx), 2 );
		    #endif
	    #elif dim == 2
            for( int k = 0; k < N; ++k ){
                int kp1 = ( k==N-1 ) ? 0 : k+1;
                int kp2 = ( k>=N-2 ) ? k-N+2 : k+2;
                int km1 = ( k==0) ? N-1: k-1;
                int km2 = ( k<2 ) ? k+N-2 : k-2;
                gradient_energy += pow( ( - f[jp2*N+k] + 8*f[jp1*N+k] - 8*f[jm1*N+k] + f[jm2*N+k] ) / (12*_dx), 2 );
                gradient_energy += pow( ( - f[j*N+kp2] + 8*f[j*N+kp1] - 8*f[j*N+km1] + f[j*N+km2] ) / (12*_dx), 2 );
            }
        #elif dim == 3
            for( int k = 0; k < N; ++k ){
                int kp1 = ( k==N-1 ) ? 0 : k+1;
                int kp2 = ( k>=N-2) ? k-N+2 : k+2;
                int km1 = ( k==0) ? N-1: k-1;
                int km2 = ( k<2 ) ? k+N-2 : k-2;
                for( int l = 0; l < N; ++l ){
                    int lp1 = ( l==N-1 ) ? 0 : l+1;
                    int lp2 = ( l>=N-2 ) ? l-N+2 : l+2;
                    int lm1 = ( l==0 ) ? N-1 : l-1;
                    int lm2 = ( l<2 ) ? l+N-2 : l-2;                        
                    gradient_energy += pow( ( - f[(jp2*N+k)*N+l] + 8*f[(jp1*N+k)*N+l] - 8*f[(jm1*N+k)*N+l] + f[(jm2*N+k)*N+l] ) / (12*_dx), 2 );
                    gradient_energy += pow( ( - f[(j*N+kp2)*N+l] + 8*f[(j*N+kp1)*N+l] - 8*f[(j*N+km1)*N+l] + f[(j*N+km2)*N+l] ) / (12*_dx), 2 );
                    gradient_energy += pow( ( - f[(j*N+k)*N+lp2] + 8*f[(j*N+k)*N+lp1] - 8*f[(j*N+k)*N+lm1] + f[(j*N+k)*N+lm2] ) / (12*_dx), 2 );
                }
            }
        #endif
    }

	return gradient_energy;
}