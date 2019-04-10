/**
 * @file parameter.hpp
 * @brief The header file written the difinitions of the parameters used in the simulation.
 * @details These words are globally declared, so beware of the confliction of the definition of local variavles.
 */
#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstddef>  // for std::size_t

/**
 * @def DIMENSION
 * @brief Simulatioin dimension.
 */
#define DIMENSION 3



/**
 * @enum Expansion
 * @brief `Expansion` parameter sets the type of cosmological expansion.
 */
enum class Expansion
{
	no_expansion,    //!< No expansion.
	rad_domination,  //!< Radiation dominated universe.
	mat_domination,  //!< Matter dominated universe.
	self_consist,    //!< Fields' energy and expansion rate are self-consistent (inflaton case).
};

extern int N;            //!< Gird size in a dimension. The total grid size is \f$ N^{dim} \f$.
extern int L;            //!< Box size of the simulation.
extern int rnd;          //!< The seed of the randum numbers.
extern int num_fields;   //!< The number of fields you use in the simulation. 
extern int num_threads;  //!< The number of threads used by OpenMP. This value  greater than 8.

extern int output_step;  //!< \f$ {\rm output\_step} \times dt \f$
extern int total_step;   //!< Total simulation time is caluculated by \f$ {\rm total_step} \times dt \f$.

extern double t0;  //!< Initial time.
extern double dt;  //!< Time step. 
extern double dx;  //!< Grid spacing calculated by \f$ L/N \f$.

extern int precision;  //!< Set the precisioin of the time evolution scheme.
extern Expansion expansion;


#endif