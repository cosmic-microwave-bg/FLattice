#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <string>
#include "calculation.hpp"


#define Logout(...) 	do { printf(__VA_ARGS__); fflush(stdout); } while(0)


void write_VTK( double* f, std::string str, int loop );

void write_status( Field* field, LeapFrog* leapfrog, Energy* energy, double** f, double t );



#endif