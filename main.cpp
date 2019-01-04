#include <chrono>
#include "parameter.hpp"
#include "utilities.hpp"



int main( int argc, char** argv )
{   
	Logout( "\n----------------------------------------------\n" );
	Logout( "            SIMULATION PARAMETERS            \n\n" );
	
	Logout( " Dimension    =  %d\n", dim );
	Logout( " Box Size     =  %d\n", L );
	Logout( " Gid Number   =  %d\n", N );
	Logout( " Initial Time =  %4.1lf\n", t0 );
	Logout( " Final Time   =  %4.2lf\n", t0 + total_step*dt );
	Logout( " dt     =  %2.4f\n", dt );
	Logout( " dx     =  %2.4f\n\n", dx );


	//--------------------------------------------------
	//       SETTING INITIAL CONDITIONS
	//--------------------------------------------------
	
	double **f, **df;
	initialize( f, df );
	write_VTK( f[0], "field", -1 );
	
	Field field( dx, num_fields, num_threads, rnd );
	LeapFrog leapfrog( precision, num_fields, num_threads, output_step, dt );

	
	//--------------------------------------------------
	//       THE TIME ITERATION LOOP
	//--------------------------------------------------

	Logout("----------------------------------------------\n");
	Logout("            STARTING COMPUTATION                \n\n");

	std::chrono::high_resolution_clock::time_point start, loop_start, current;
	std::chrono::milliseconds elapsed;

	start = std::chrono::high_resolution_clock::now();
    
    for( int loop = 0; loop < max_loop; ++loop )
	{	
	    double t = t0 + loop*output_step*dt;
	    loop_start = std::chrono::high_resolution_clock::now();

		//leapfrog.evolution( &field, f, df );
		leapfrog.evolution_expansion( &field, f, df, t );

		Energy energy( &field, &leapfrog, f, df, num_fields, dx );
		write_VTK( f[0], "field", loop );
		//write_VTK( energy.value[0], "energy", loop );
		write_status( &field, &leapfrog, &energy, f, t );

        current = std::chrono::high_resolution_clock::now();
	    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( current - loop_start );
        Logout( " Loop %d/%d: %2.3f s \n", loop+1, max_loop, elapsed.count()*1.e-3 );
    }

    current = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( current - start );
    Logout( " Total time: %2.3f s \n", elapsed.count()*1.e-3 );
    
	finalize( f, df );

	Logout( "\n PROGRAMM FINISHED \n" );
	Logout( "----------------------------------------------\n\n" );
}