#include "parameter.hpp"


int N = 256;
int L = 10;
int rnd = 1;
int num_fields  = 2;
int num_threads = 2;

int output_step = 1e+1;
int total_step  = 1e+3;
int max_loop    = total_step/output_step;

double t0 = 1.;
double dt = 1.e-2;
double dx = 1.* L/N;

int expansion = 0; // 0: no expanxion, 1: self-consistent, 2: radiation dominant, 3: matter dominant
int precision = 2;

bool restart = false;  // cannot use yet

