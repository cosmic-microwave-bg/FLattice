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

int output_step;
int total_step;

double t0;
double dt;
double dx;

int precision;
Expansion expansion;

double K;
double meff2;

double a;
double da;
double dda;
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

    output_step = j["output_step"];
    total_step  = j["total_step"];

    t0 = j["t0"];
    dt = j["dt"];
    dx = 1.* L/N;

    precision = j["precision"];
    expansion = j["expansion"].get<Expansion>();

    K = j["K"];
    meff2 = j["meff2"];

    a   = 1;
    da  = 0;
    dda = 0;
    t   = t0;
}
