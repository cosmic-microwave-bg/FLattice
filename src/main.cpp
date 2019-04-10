#include <random>
#include <memory>
#include "parameter.hpp"
#include "utilities.hpp"
#include "stopwatch.hpp"


void initialize( double** f, double** df )
{   
	std::mt19937 mt( rnd );
	std::uniform_real_distribution<> rand( -1.e-1, 1.e-1 );

    for( int i = 0; i < num_fields; ++i ){
        #pragma omp parallel for schedule( static ) num_threads ( num_threads )
        for( int j = 0; j < N; ++j ){	
            for( int k = 0; k < N; ++k ){
				int idx = j*N + k;
				f[i][idx] = 10*( 1 + rand(mt) );
				df[i][idx] = 0;
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

	auto model    = std::make_shared<Model>();
	auto leapfrog = std::make_shared<LeapFrog>(model);

	initialize( f, df );

	//write_VTK( f[0], "model", -1 );
	CalculateEnergy energy( model, f, df );
	// write_VTK( energy.value[0], "energy", -1 );
	//write_status( &model, &leapfrog, &energy, f, t0 );

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
		
		CalculateEnergy energy( model, f, df );
		// write_VTK( f[0], "model", loop );
		// write_VTK( energy.value[0], "energy", loop );
		//write_status( &model, &leapfrog, &energy, f, t+output_step*dt );

        Logout( " Loop %d/%d: %2.3f s \n", loop, max_loop, stopwatch.lap() );
    }

    Logout( " Total time: %2.3f s \n", stopwatch.end() );
    
	deleteData<double>( f);
	deleteData<double>(df);

	Logout( "\n PROGRAMM FINISHED \n" );
	Logout( "----------------------------------------------\n\n" );
}