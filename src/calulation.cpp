#include <cmath>

#include "calculation.hpp"
#include "parameter.hpp"



CalculateEnergy::CalculateEnergy( std::shared_ptr<Model> model, double** f, double** df ): CalculateBase()
{
    double a  = 1;//leapfrog->a ();
    double da = 1;//leapfrog->da();
    
    for( int i = 0; i < num_fields; ++i ){
        // #pragma omp parallel for schedule( static ) num_threads ( num_threads )
        for( int j = 0; j < N; ++j ){
            if( expansion == Expansion::no_expansion )
            {
                switch( DIMENSION ){
                    case 1:
                        int idx = j;
                        _data[i][idx] = pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,4)) + model->aV(f, a, i ,idx);
                        _average[i] += _data[i][idx];
                        break;
                    case 2:
                        for( int k = 0; k < N; ++k ){
                            int idx = j*N + k;
                            _data[i][idx] = pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,4)) + model->aV(f, a, i ,idx);
                            _average[i] += _data[i][idx];
                        }
                        break;
                    case 3:
                        for( int k = 0; k < N; ++k ){
                            for( int l = 0; l < N; ++l ){
                                int idx = (j*N + k)*N + l;
                                _data[i][idx] = pow(df[i][idx]*a - f[i][idx]*da, 2)/(2*pow(a,4)) + model->aV(f, a, i ,idx);
                                _average[i] += _data[i][idx];
                            }
                        }
                        break;
                }
            }
            else
            {
                switch( DIMENSION ){
                    case 1:
                        int idx = j;
                        _data[i][idx] = pow(df[i][idx], 2)/2 + model->V(f, i ,idx);
                        _average[i] += _data[i][idx];
                        break;
                    case 2:
                        for( int k = 0; k < N; ++k ){
                            int idx = j*N + k;
                            _data[i][idx] = pow(df[i][idx], 2)/2 + model->V(f, i ,idx);
                            _average[i] += _data[i][idx];
                        }
                        break;
                    case 3:
                        for( int k = 0; k < N; ++k ){
                            for( int l = 0; l < N; ++l ){
                                int idx = (j*N + k)*N + l;
                                _data[i][idx] = pow(df[i][idx], 2)/2 + model->V(f, i ,idx);
                                _average[i] += _data[i][idx];
                            }
                        }
                        break;
                }
            }
        }
    }
    
    for( int n = 0; n < num_fields; ++n ){
        if( expansion == Expansion::no_expansion ) _average[n] += gradient_energy(f[n]);
        else                                       _average[n] += gradient_energy(f[n])/ pow(a,4);
        for( int j = 0; j < DIMENSION; ++j ) _average[n] /= N;
    }
}


double CalculateEnergy::gradient_energy( double* f )
{
    double gradient_energy = 0;
    
    #pragma omp parallel for reduction(+:gradient_energy) schedule( static ) num_threads ( num_threads )
    for( int j = 0; j < N; ++j ){
        int jp1 = (j == N-1)?     0: j+1;
        int jp2 = (j >= N-2)? j-N+2: j+2;
        #ifdef SPHERICAL_SYM
        int jm1 = abs( j-1 );
        int jm2 = abs( j-2 );
        #else
        int jm1 = (j == 0)?   N-1: j-1;
        int jm2 = (j <  2)? j+N-2: j-2;
        #endif
        #if DIMENSION == 1
		    #ifdef SPHERICAL_SYM
            int idx = j;
            if( idx == 0 ){
                grad[0] = 0;
            }else if( idx == N-2 ){
                grad[0] = ( f[jp1] - f[jm1] ) / (2*dx);
            }else if( idx == N-1 ){
                grad[0] = ( f[jm2] - 4*f[jm2] + 3*f[idx]  ) / (2*dx);
            }else{
                grad[0] = ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx);
            }
            #else
            gradient_energy +=  pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx), 2 );
		    #endif
	    #elif DIMENSION == 2
            for( int k = 0; k < N; ++k ){
                int kp1 = (k == N-1)?     0: k+1;
                int kp2 = (k >= N-2)? k-N+2: k+2;
                int km1 = (k ==   0)?   N-1: k-1;
                int km2 = (k <    2)? k+N-2: k-2;
                gradient_energy += pow( (- f[jp2*N+k] + 8*f[jp1*N+k] - 8*f[jm1*N+k] + f[jm2*N+k]) / (12*dx), 2 );
                gradient_energy += pow( (- f[j*N+kp2] + 8*f[j*N+kp1] - 8*f[j*N+km1] + f[j*N+km2]) / (12*dx), 2 );
            }
        #elif DIMENSION == 3
            for( int k = 0; k < N; ++k ){
                int kp1 = (k == N-1)?     0: k+1;
                int kp2 = (k >= N-2)? k-N+2: k+2;
                int km1 = (k ==   0)?   N-1: k-1;
                int km2 = (k <    2)? k+N-2: k-2;
                for( int l = 0; l < N; ++l ){
                    int lp1 = (l == N-1)?     0: l+1;
                    int lp2 = (l >= N-2)? l-N+2: l+2;
                    int lm1 = (l ==   0)?   N-1: l-1;
                    int lm2 = (l <    2)? l+N-2: l-2;                        
                    gradient_energy += pow( (- f[(jp2*N+k)*N+l] + 8*f[(jp1*N+k)*N+l] - 8*f[(jm1*N+k)*N+l] + f[(jm2*N+k)*N+l]) / (12*dx), 2 );
                    gradient_energy += pow( (- f[(j*N+kp2)*N+l] + 8*f[(j*N+kp1)*N+l] - 8*f[(j*N+km1)*N+l] + f[(j*N+km2)*N+l]) / (12*dx), 2 );
                    gradient_energy += pow( (- f[(j*N+k)*N+lp2] + 8*f[(j*N+k)*N+lp1] - 8*f[(j*N+k)*N+lm1] + f[(j*N+k)*N+lm2]) / (12*dx), 2 );
                }
            }
        #endif
    }

	return gradient_energy / 2;
}