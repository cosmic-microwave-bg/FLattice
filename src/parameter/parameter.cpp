/**
 * @file parameter.cpp
 * @brief See the explanation of the coresponding header file.
 */
#include "parameter.hpp"


int N;
int L;
int rnd;
int num_fields;
int num_threads;

double dt;
int output_step;
int total_step;

int A;
int D;

Expansion expansion;
double c;
double t0;
double R;

int precision;

double K;
double meff2;

double dx;
double a, da;
double t;


void setParameters()
{
    nlohmann::json j;
    std::ifstream input("../src/parameter/parameter.json");
    input >> j;

    N = j["N"];
    L = j["L"];
    rnd = j["rnd"];
    num_fields  = j["num_fields"];
    num_threads = j["num_threads"];

    dt = j["dt"];
    output_step = j["output_step"];
    total_step  = j["total_step"];

    A = j["A"];
    D = j["D"];

    expansion = j["expansion"].get<Expansion>();
    c = j["c"];
    t0 = j["t0"];
    R = j["F/Mp"];

    precision = j["precision"];

    K = j["K"];
    meff2 = j["meff2"];

    dx = 1.* L/N;
    a   = 1;
    da  = 0;
    t   = t0;
}
