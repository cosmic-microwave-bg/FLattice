/**
 * @file    parameter.hpp
 * @brief   The header file written the difinitions of the parameters used in the simulation.
 * @details These words are globally declared, so beware of the confliction of the definition of local variavles.
 *          This file will be changed to .json style in the future update
 */
#ifndef _PARAMETER_H_
#define _PARAMETER_H_


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstddef>  // for std::size_t
#include "nlohmann/json.hpp"

/**
 * @def DIMENSION
 * @brief Simulation dimension.
 */
#define DIMENSION 3
#define Logout(...)    do { printf(__VA_ARGS__); fflush(stdout); } while(0)



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
NLOHMANN_JSON_SERIALIZE_ENUM( Expansion, {
    {Expansion::no_expansion, "no expansion"},
    {Expansion::rad_domination, "radiation dominated universe"},
    {Expansion::mat_domination, "matter dominated universe"},
    {Expansion::self_consist, "self consistent universe"},
})

extern int A;
extern int D;
extern Expansion expansion;
extern double R;  //!< Field normalization to reduced Planck mass
extern double meff2;


extern int N;            //!< Gird size in a dimension. The total grid size is \f$ N^{dim} \f$.
extern int L;            //!< Box size of the simulation.
extern int num_fields;   //!< The number of fields you use in the simulation. 
extern int num_threads;  //!< The number of threads used by OpenMP. This value  greater than 8.

extern double dt;  //!< Time step. 
extern int output_step;  //!< \f$ {\rm output\_step} \times dt \f$
extern int total_step;   //!< Total simulation time is caluculated by \f$ {\rm total_step} \times dt \f$.

extern int precision;  //!< Set the precisioin of the time evolution scheme.


extern int rnd;          //!< The seed of the randum numbers.
extern double c;
extern double t0;  //!< Initial time.
extern std::vector<double> ini_amp;
extern std::vector<double> ini_vel;


extern double lambda;


extern double dx;  //!< Grid spacing calculated by \f$ L/N \f$.
extern double a, da;
extern double t;


void setParameters();



#endif