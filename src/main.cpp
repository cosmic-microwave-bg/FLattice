/**
 * @file  main.cpp
 * @brief Main file of Flattice.
 */
#include <random>
#include <memory>
#include <fftw3.h>
#include "status.hpp"
#include "stopwatch.hpp"



inline void DFT_c2r( double** f )
{
	fftw_plan p;
	double* out;
    fftw_complex *in;
	size_t in_size;
	  
	switch( DIMENSION )
	{
		case 1:
			in_size = sizeof(fftw_complex) * (N/2+1);
			in  = (fftw_complex*)fftw_malloc( in_size );
			out = new double [N]();
			p = fftw_plan_dft_c2r_1d( N, in, out, FFTW_ESTIMATE );
			break;
		case 2:
			in_size = sizeof(fftw_complex) * N * (N/2+1);
			in  = (fftw_complex*)fftw_malloc( in_size );
			out = new double [N*N]();
			p = fftw_plan_dft_c2r_2d( N, N, in, out, FFTW_ESTIMATE );
			break;
		case 3:
			in_size = sizeof(fftw_complex) * N * N * (N/2+1);
			in  = (fftw_complex*)fftw_malloc( in_size );
			out = new double [N*N*N]();
			p = fftw_plan_dft_c2r_3d( N, N, N, in, out, FFTW_ESTIMATE );
			break;
	}

	for( int n = 0; n < num_fields; ++n )
	{
		std::mt19937 mt( rnd );
		std::uniform_real_distribution<> rand( 0, 2*M_PI );
		
		// Create input data
		switch( DIMENSION ){
			case 1:
				#pragma omp parallel for schedule( static ) num_threads( num_threads )
				for( int i = 0; i < N/2+1; ++i ){
					int idx = i;
					double phase = rand(mt);
					in[idx][0] = f[n][idx] * cos( phase );
					in[idx][1] = f[n][idx] * sin( phase );
				}
				break;
			case 2:
				#pragma omp parallel for schedule( static ) num_threads( num_threads )
				for( int i = 0; i < N; ++i ){
					for( int j = 0; j < N/2+1; ++j ){
							int idx = i*N + j;
							double phase = rand(mt);
							in[idx][0] = f[n][idx] * cos( phase );
							in[idx][1] = f[n][idx] * sin( phase );
					}
				}
				break;
			case 3:
				#pragma omp parallel for schedule( static ) num_threads( num_threads )
				for( int j = 0; j < N; ++j ){
					for( int k = 0; k < N; ++k ){
							for( int l = 0; l < N/2+1; ++l ){
									int idx = (j*N + k)*N + l;
									double phase = rand(mt);
									in[idx][0] = f[n][idx] * cos( phase );
									in[idx][1] = f[n][idx] * sin( phase );
							}
					}
				}
				break;
		}
	
		fftw_execute(p);

		// Set output data
		#pragma omp parallel for schedule( static ) num_threads( num_threads )
		for( int j = 0; j < N; ++j ){
			switch( DIMENSION ){
				case 1:
					int idx = j;
					f[n][idx] = out[idx]/N;
					break;
				case 2:
					for( int k = 0; k < N; ++k ){
						int idx = j*N + k;
						f[n][idx] = out[idx]/(N*N);
					}
				break;
				case 3:
					for( int k = 0; k < N; ++k ){
						for( int l = 0; l < N; ++l ){
							int idx = (j*N + k)*N + l;
							f[n][idx] = out[idx]/(N*N*N);
						}
					}
				break;
			}
		}
	}

	if( p ) fftw_destroy_plan(p);
	if( in ) fftw_free(in);
	delete[] out;
}


void initialize( double** f, double** df )
{   
	std::mt19937 mt( rnd );
	std::uniform_real_distribution<> rand( -1.e-2, 1.e-2 );

    for( int n = 0; n < num_fields; ++n ){
        #pragma omp parallel for schedule( static ) num_threads ( num_threads )
        for( int i = 0; i < N; ++i ){	
            for( int j = 0; j < N; ++j ){
            	for( int k = 0; k < N; ++k ){
					int idx = (i*N+j)*N+k;
					f[n][idx]  = 10*( 1 + rand(mt) );
					df[n][idx] = 0;
				}
            }
        }
    }

	//DFT_c2r( f );

}



int main( int argc, char** argv )
{   
	Logout( "\n----------------------------------------------\n" );
	Logout( "            SIMULATION PARAMETERS            \n\n" );
	
	Logout( " Dimension    =  %d\n",  DIMENSION );
	Logout( " Box Size     =  %d\n",  L );
	Logout( " Gid Number   =  %d\n",  N );
	Logout( " Initial Time =  %4.1lf\n", t0 );
	Logout( " Final Time   =  %4.2e\n",  t0 + total_step*dt );
	Logout( " dt     =  %2.2e\n", dt );
	Logout( " dx     =  %2.2e\n", dx );
	Logout( " Number of fields   =  %d\n", num_fields );
	Logout( " Number of threads  =  %d\n\n", num_threads );


	//-------------------------------------------------
	//       SETTING INITIAL CONDITIONS
	//-------------------------------------------------
	
	double **f  = allocateData<double>( num_fields, N, DIMENSION );
	double **df = allocateData<double>( num_fields, N, DIMENSION );

	initialize( f, df );

	auto model    = std::make_shared<Model>();
	auto leapfrog = std::make_shared<LeapFrog>(model);
	auto energy   = std::make_shared<CalculateEnergy>(model, f, df);
	auto status   = std::make_shared<Status>(f, energy);

	energy->update();
	//energy.write()

	status->update();
	status->write();
	// write_VTK( energy.value[0], "energy", -1 );
	// write_VTK( f[0], "model", -1 );

	//-------------------------------------------------
	//           TIME ITERATION LOOP
	//-------------------------------------------------

	Logout("----------------------------------------------\n");
	Logout("            STARTING COMPUTATION                \n\n");

	Stopwatch stopwatch;
	int max_loop = total_step/output_step;

    for( int loop = 1; loop <= max_loop; ++loop )
	{	
		leapfrog->evolution( f, df );
		
		energy->update(loop);
		status->update(loop);
		status->write (loop);
		// write_VTK( f[0], "model", loop );
		// write_VTK( energy.value[0], "energy", loop );

        Logout( " Loop %d/%d: %2.3f s \n", loop, max_loop, stopwatch.lap() );
    }

    Logout( " Total time: %2.3f s \n", stopwatch.end() );
    
	deleteData<double>( f);
	deleteData<double>(df);

	Logout( "\n PROGRAMM FINISHED \n" );
	Logout( "----------------------------------------------\n\n" );
}