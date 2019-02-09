#include <cmath>
#include <random>
#include "field.hpp"
#include "utilities.hpp"


void initialize( double**& f, double**& df )
{   
	f = new double* [num_fields];
	df = new double* [num_fields];

	switch( dim )
	{
		case 1:
			f[0] = new double [num_fields*N];
			df[0] = new double [num_fields*N];
			for( int i = 0; i < num_fields; ++i )
			{
				f[i] = f[0] + i*N;
				df[i] = df[0] + i*N;
			}
			break;
		case 2:
			f[0] = new double [num_fields*N*N];
			df[0] = new double [num_fields*N*N];
			for( int i = 0; i < num_fields; ++i )
			{
				f[i] = f[0] + i*N*N;
				df[i] = df[0] + i*N*N;
			}
			break;
		case 3:
			f[0] = new double [num_fields*N*N*N];
			df[0] = new double [num_fields*N*N*N];
			for( int i = 0; i < num_fields; ++i )
			{
				f[i] = f[0] + i*N*N*N;
				df[i] = df[0] + i*N*N*N;
			}
			break;
		default:
			std::cout << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
			exit(1);
	}

	std::mt19937 mt( rnd );
	std::uniform_real_distribution<> rand( -1.e-2, 1.e-2 );

    for( int i = 0; i < num_fields; ++i ){
        #pragma omp parallel for schedule( static ) num_threads ( num_threads )
        for( int j = 0; j < N; ++j ){	
            for( int k = 0; k < N; ++k ){
				int idx = j*N + k;
				if( idx = 10 ) f[i][idx] = 10;
				else f[i][idx] = 0;
				df[i][idx] = 0;
            }
        }
    }

	DFT_c2r( f );

}

void finalize( double**& f, double**& df )
{
	delete [] f[0];
	delete [] df[0];
    delete [] f;
	delete [] df;
}


Field::Field( double dx, int num_fields, int num_threads, int rnd ):
                _dx(dx), _num_fields(num_fields), _num_threads(num_threads), _rnd(rnd)
{
    _average = new double [_num_fields];
	_variance = new double [_num_fields];
}


double Field::laplacian( double* f, int j, int k, int l ) const
{	
    int jp1 = ( j == N-1 ) ? 0 : j+1;
	int jp2 = ( j >= N-2 ) ? j-N+2 : j+2;
	#ifdef SPHERICAL_SYM
	if( dim != 1 ){
		Logout( "Error: Dimension must be 1 when you define SPHERICAL_SYM with ABC. \n" );
		exit(1);
	}
	int jm1 = abs( j-1 );
	int jm2 = abs( j-2 );
	#else
	int jm1 = ( j == 0 ) ? N-1 : j-1;
	int jm2 = ( j < 2 ) ? j+N-2 : j-2;
	#endif
	
	#if dim == 1
	    int idx = j;
		#ifdef SPHERICAL_SYM
            if( idx == 0 ){ // r = 0では1/r*df/dr = 0 より別に処理
                return (- f[jp2] + 16*f[jp1] - 30*f[idx] + 16*f[jm1] - f[jm2]) / (12*_dx*_dx);
            }else if( idx == N-2 ){ //境界1つ手前は中央差分2次で計算
                return (f[jp1] - 2*f[idx] + f[jm1]) / (_dx*_dx) + 2*gradient(f, 0, i, idx, 0, 0) / (idx*_dx);
            }else if( idx == N-1 ){ //境界では後退差分2次で計算、時間発展はABC
                return (df[jm2] - 4*df[jm1] + 3*df[idx]) / (2*_dx) + df[idx] / (idx*_dx);
            }else{
                return (- f[jp2] + 16*f[jp1] - 30*f[idx] + 16*f[jm1] - f[jm2]) / (12*_dx*_dx) + 2*gradient(f, 0, i, idx, 0, 0) / (idx*_dx);
            }
		#else
		    return (- f[jp2] + 16*f[jp1] - 30*f[idx] + 16*f[jm1] - f[jm2]) / (12*_dx*_dx);
		#endif
	#elif dim == 2
		int kp1 = ( k == N-1 ) ? 0 : k+1;
		int kp2 = ( k >= N-2 ) ? k-N+2 : k+2;
		int km1 = ( k == 0 ) ? N-1: k-1;
		int km2 = ( k < 2 ) ? k+N-2 : k-2;

        int idx = j*N + k;
        return ( (- f[jp2*N+k] + 16*f[jp1*N+k] - 30*f[idx] + 16*f[jm1*N+k] - f[jm2*N+k]) 
                + (- f[j*N+kp2] + 16*f[j*N+kp1] - 30*f[idx] + 16*f[j*N+km1] - f[j*N+km2]) ) / (12*_dx*_dx);
	#elif dim == 3
		int kp1 = ( k == N-1 ) ? 0 : k+1;
		int kp2 = ( k >= N-2 ) ? k-N+2 : k+2;
		int km1 = ( k == 0 ) ? N-1: k-1;
		int km2 = ( k < 2 ) ? k+N-2 : k-2;
		
		int lp1 = ( l == N-1 ) ? 0 : l+1;
		int lp2 = ( l >= N-2 ) ? l-N+2 : l+2;
		int lm1 = ( l == 0 ) ? N-1 : l-1;
		int lm2 = ( l < 2 ) ? l+N-2 : l-2;

        int idx = (j*N+k)*N + l;
        return ( (- f[(jp2*N+k)*N+l] + 16*f[(jp1*N+k)*N+l] - 30*f[idx] + 16*f[(jm1*N+k)*N+l] - f[(jm2*N+k)*N+l])
                + (- f[(j*N+kp2)*N+l] + 16*f[(j*N+kp1)*N+l] - 30*f[idx] + 16*f[(j*N+km1)*N+l] - f[(j*N+km2)*N+l])
                + (- f[(j*N+k)*N+lp2] + 16*f[(j*N+k)*N+lp1] - 30*f[idx] + 16*f[(j*N+k)*N+lm1] - f[(j*N+k)*N+lm2]) ) / (12*_dx*_dx);
	#endif
}


double Field::gradient( double* f, int d, int j, int k, int l ) const
{
    double grad[dim];
	
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
		}else if( idx == 1 || idx == N-2 ){
			grad[0] = ( f[jp1] - f[jm1] ) / (2*_dx);
		}else if( idx == N-1 ){
			grad[0] = ( f[jm2] - 4*f[jm2] + 3*f[idx]  ) / (2*_dx);
		}else{
			grad[0] = ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*_dx);
		}
		#else
		grad[0] = ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*_dx);
		#endif
	#elif dim == 2
		int kp1 = ( k==N-1 ) ? 0 : k+1;
		int kp2 = ( k>=N-2) ? k-N+2 : k+2;
		int km1 = ( k==0) ? N-1: k-1;
		int km2 = ( k<2 ) ? k+N-2 : k-2;
		grad[0] = ( - f[jp2*N+k] + 8*f[jp1*N+k] - 8*f[jm1*N+k] + f[jm2*N+k] ) / (12*_dx);
		grad[1] = ( - f[j*N+kp2] + 8*f[j*N+kp1]  - 8*f[j*N+km1] + f[j*N+km2] ) / (12*_dx);
	#elif dim == 3
		int kp1 = ( k==N-1 ) ? 0 : k+1;
		int kp2 = ( k>=N-2 ) ? k-N+2 : k+2;
		int km1 = ( k==0 ) ? N-1: k-1;
		int km2 = ( k<2 ) ? k+N-2 : k-2;
		
		int lp1 = ( l==N-1 ) ? 0 : l+1;
		int lp2 = ( l>=N-2 ) ? l-N+2 : l+2;
		int lm1 = ( l==0 ) ? N-1 : l-1;
		int lm2 = ( l<2 ) ? l+N-2 : l-2;
		grad[0] = ( - f[(jp2*N+k)*N+l] + 8*f[(jp1*N+k)*N+l] - 8*f[(jm1*N+k)*N+l] + f[(jm2*N+k)*N+l] ) / (12*_dx);
		grad[1] = ( - f[(j*N+kp2)*N+l] + 8*f[(j*N+kp1)*N+l] - 8*f[(j*N+km1)*N+l] + f[(j*N+km2)*N+l] ) / (12*_dx);
		grad[2] = ( - f[(j*N+k)*N+lp2] + 8*f[(j*N+k)*N+lp1] - 8*f[(j*N+k)*N+lm1] + f[(j*N+k)*N+lm2] ) / (12*_dx);
	#endif
	
	return grad[d];
}


double Field::gradient_energy( double* f ) const
{
    double gradient_energy = 0;
    
    #pragma omp parallel for schedule( static ) num_threads ( num_threads )
    for( int j = 0; j < N; ++j ){
        int jp1 = ( j ==N-1) ? 0 : j+1;
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
                int kp2 = ( k>=N-2) ? k-N+2 : k+2;
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

double Field::potential_energy( double* f, double a ) const
{
	double potential_energy = 0;
    
    #pragma omp parallel for schedule( static ) num_threads ( num_threads )
	for( int j = 0; j < N; ++j ){
        #if dim == 1
            potential_energy +=  pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*_dx), 2 );
	    #elif dim == 2
            for( int k = 0; k < N; ++k )
			{
				int idx = j*N + k;
                //potential_energy += aV( f, a, 0, idx );
            }
        #elif dim == 3
            for( int k = 0; k < N; ++k ){
                for( int l = 0; l < N; ++l ){
					int idx = ( j*N + k)*N + l;
					potential_energy += f[idx]*f[idx]/(2*a*a);
				}
            }
        #endif
	}
	return potential_energy;
}


double Field::average( double* f, int i )
{
	_average[i] = 0;

	for( int j = 0; j < N; ++j ){
		switch( dim )
		{
			case 1:
				int idx = j;
				_average[i] += f[idx];
				break;
			case 2:
				for( int k = 0; k < N; ++k ){
					int idx = j*N + k;
					_average[i] += f[idx];
				}
				break;
			case 3:
				for( int k = 0; k < N; ++k ){
					for( int l = 0; l < N; ++l ){
						int idx = (j*N + k)*N + l;
						_average[i] += f[idx];
					}
				}
				break;
		}
	}
    for( int j = 0; j < dim; ++j ) _average[i] /= N;
	
	return _average[i];
}


double Field::variance( double* f, int i )
{
	_variance[i] = 0;

	for( int j = 0; j < N; ++j ){
		switch( dim ){
			case 1:
				int idx = j;
				_variance[i] += pow( f[idx] - _average[i], 2 );
				break;
			case 2:
				for( int k = 0; k < N; ++k ){
					int idx = j*N + k;
					_variance[i] += pow( f[idx] - _average[i], 2 );
				}
				break;
			case 3:
				for( int k = 0; k < N; ++k ){
					for( int l = 0; l < N; ++l ){
						int idx = (j*N + k)*N + l;
						_variance[i] += pow( f[idx] - _average[i], 2 );
					}
				}
				break;
		}
	}
    for( int j = 0; j < dim; ++j ) _variance[i] /= N;
	
	return _variance[i];
}