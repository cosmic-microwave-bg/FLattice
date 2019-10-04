/**
 * @file  main.cpp
 * @brief Main file of Flattice.
 */
#include <random>
#include <memory>
#include "simulator.hpp"
#include "stopwatch.hpp"



int main()
{   
    setParameters();
    Simulator simulator(num_fields, N, DIMENSION);

    //-------------------------------------------------
    //       SETTING INITIAL CONDITIONS
    //-------------------------------------------------
    
    simulator.initializeField(rnd);
    simulator.setModel("harmonic_oscillator");
    simulator.setEvolutionScheme("leapfrog", precision, expansion);
    
    simulator.addPhysicalQuantites("energy");
    simulator.addPhysicalQuantites("charge");
    simulator.calculatePhysicalQuantities();

    simulator.writeFields();
    simulator.writePhysicalQuantities();
    simulator.writeStatus();

    //-------------------------------------------------
    //           TIME ITERATION LOOP
    //-------------------------------------------------

    Logout("----------------------------------------------\n");
    Logout("            STARTING COMPUTATION                \n\n");

    Stopwatch stopwatch;
    int max_loop = total_step/output_step;

    for( int loop = 1; loop <= max_loop; ++loop )
    {
        simulator.run(output_step);
        simulator.calculatePhysicalQuantities();

        t = t0 + loop*(dt*output_step);

        simulator.writeFields(loop);
        simulator.writePhysicalQuantities(loop);
        simulator.writeStatus();

        Logout( " Loop %d/%d: %2.3f s \n", loop, max_loop, stopwatch.lap() );
    }
    Logout( " Total time: %2.3f s \n", stopwatch.end() );
    

    Logout( "\n PROGRAMM FINISHED \n" );
    Logout( "----------------------------------------------\n\n" );
}