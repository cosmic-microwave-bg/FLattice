#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define dim 3


extern int N;
extern int L;
extern int rnd;
extern int num_fields;
extern int num_threads;

extern int output_step;
extern int total_step;
extern int max_loop;

extern double t0;
extern double dt;
extern double dx;

extern int expansion;
extern int precision;
extern bool restart;



#endif