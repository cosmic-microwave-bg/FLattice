/**
 * @file parameter.cpp
 * @brief See the explanation of the coresponding header file.
 */
#include "parameter.hpp"


int N = 64;
int L = 10;
int rnd = 1;
int num_fields  = 2;
int num_threads = 2;

int output_step = 1e+2;
int total_step  = 5e+2;

double t0 = 1.;
double dt = 1.e-2;
double dx = 1.* L/N;

int precision = 4;
Expansion expansion = Expansion::rad_domination;