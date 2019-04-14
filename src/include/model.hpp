/**
 * @file model.hpp
 * @brief class `Model` is written where model parameters and the potential is defined.
 * 
 * Only you have to do is to change the potential and derivative of potential of model class.
 */
#ifndef _MODEL_H_
#define _MODEL_H_

#include <cmath>
#include "parameter.hpp"



template <typename dataType>
dataType** allocateData( int num_fields, std::size_t N, int dimension )
{
	dataType** data = new dataType* [num_fields];
	switch( dimension )
	{
		case 1:
			data[0] = new dataType [num_fields*N]();
			for( int n = 0; n < num_fields; ++n ) data[n] = data[0] + n*N;
			break;
		case 2:
			data[0] = new dataType [num_fields*N*N]();
			for( int n = 0; n < num_fields; ++n ) data[n] = data[0] + n*N*N;
			break;
		case 3:
			data[0] = new dataType [num_fields*N*N*N]();
			for( int n = 0; n < num_fields; ++n ) data[n] = data[0] + n*N*N*N;
			break;
		default:
			std::cerr << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
			exit(1);
	}
	return data;
}


template <typename dataType>
void deleteData( dataType**& data )
{
	delete [] data[0];
	delete [] data;
}


/**
 * @class Model
 * @brief You set the model potential (and parameters) in this class model.
 */
class Model
{
	public:
		/**
		 * @brief   Potential
		 * @details When you simulate in the Minkowski spacetime, just write the potential \f$ V(\phi) \f$.
		 *          When you simulate in the FRW spacetime, write the rescaled potential \f$ V(\bar{\phi}) \f$
		 *          where \f$ \bar{\phi} = \phi / a \f$
		 */
		double V   ( double** f, int n, int idx, double a=1 ) const { return pow(f[n][idx], 2)/(2*a*a); }
		double dV  ( double** f, int n, int idx, double a=1 ) const { return     f[n][idx]    /a; }

		double gradientEnergy  ( double* f ) const
		{
			double gradient_energy = 0;
			
			#pragma omp parallel for reduction (+:gradient_energy) schedule( static ) num_threads ( num_threads )
			for( int j = 0; j < N; ++j ){
				int jp1 = (j == N-1)?     0: j+1;
				int jp2 = (j >= N-2)? j-N+2: j+2;
				int jm1 = (j ==   0)?   N-1: j-1;
				int jm2 = (j <    2)? j+N-2: j-2;
				#if   DIMENSION == 1
					gradient_energy +=  pow( ( - f[jp2]  + 8*f[jp1]  - 8*f[jm1] + f[jm2] ) / (12*dx), 2 );
				#elif DIMENSION == 2
					for( int k = 0; k < N; ++k ){
						int kp1 = (k == N-1)?     0: k+1;
						int kp2 = (k >= N-2)? k-N+2: k+2;
						int km1 = (k ==   0)?   N-1: k-1;
						int km2 = (k <    2)? k+N-2: k-2;
						gradient_energy += pow( ( - f[jp2*N+k] + 8*f[jp1*N+k] - 8*f[jm1*N+k] + f[jm2*N+k] ) / (12*dx), 2 );
						gradient_energy += pow( ( - f[j*N+kp2] + 8*f[j*N+kp1] - 8*f[j*N+km1] + f[j*N+km2] ) / (12*dx), 2 );
					}
				#elif DIMENSION == 3
					for( int k = 0; k < N; ++k ){
						int kp1 = (k == N-1)?     0: k+1;
						int kp2 = (k >= N-2)? k-N+2: k+2;
						int km1 = (k ==   0)?   N-1: k-1;
						int km2 = (k <    2)? k+N-2: k-2;
						for( int l = 0; l < N; ++l ){
							int lp1 = (l == N-1) ?     0: l+1;
							int lp2 = (l >= N-2) ? l-N+2: l+2;
							int lm1 = (l ==   0) ?   N-1: l-1;
							int lm2 = (l <    2) ? l+N-2: l-2;                        
							gradient_energy += pow( ( - f[(jp2*N+k)*N+l] + 8*f[(jp1*N+k)*N+l] - 8*f[(jm1*N+k)*N+l] + f[(jm2*N+k)*N+l] ) / (12*dx), 2 );
							gradient_energy += pow( ( - f[(j*N+kp2)*N+l] + 8*f[(j*N+kp1)*N+l] - 8*f[(j*N+km1)*N+l] + f[(j*N+km2)*N+l] ) / (12*dx), 2 );
							gradient_energy += pow( ( - f[(j*N+k)*N+lp2] + 8*f[(j*N+k)*N+lp1] - 8*f[(j*N+k)*N+lm1] + f[(j*N+k)*N+lm2] ) / (12*dx), 2 );
						}
					}
				#endif
			}
			for( int d = 0; d < DIMENSION; ++d ) gradient_energy /= N;
			return gradient_energy;
		}

		double potential_energy ( double** f, int n, double a ) const
		{
			double potential_energy = 0;

			#pragma omp parallel for reduction(+:potential_energy) schedule( static ) num_threads ( num_threads )
			for( int i = 0; i < N; ++i ){
				#if   DIMENSION == 1
					int idx = i;
					potential_energy += V(f, n, idx, a);
				#elif DIMENSION == 2
					for( int j = 0; j < N; ++j ){
						int idx = i*N+j;
						potential_energy += V(f, n, idx, a);
					}
				#elif DIMENSION == 3
					for( int j = 0; j < N; ++j )
						for( int k = 0; k < N; ++k ){
							int idx = (i*N+j)*N+k;
							potential_energy += V(f, n, idx, a);
						}
				#endif
			}

			for( int d = 0; d < DIMENSION; ++d ) potential_energy /= N;
			return potential_energy;
		}
};



#endif