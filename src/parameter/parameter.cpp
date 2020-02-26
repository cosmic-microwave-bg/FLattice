/**
 * @file parameter.cpp
 * @brief See the explanation of the coresponding header file.
 */
#include <typeinfo>
#include "parameter.hpp"


int A;
int D;
Expansion expansion;
double R;
double meff2;


int N;
int L;
int num_fields;
int num_threads;

double dt;
int output_step;
int total_step;

int precision;


int rnd;
double c;
double t0;
std::vector<double> ini_amp;
std::vector<double> ini_vel;


double lambda;


double dx;
double a, da;
double t;


void setParameters()
{
    nlohmann::json j;
    std::ifstream input("../src/parameter/parameter.json");
    input >> j;


    A = j["A"];
    D = j["D"];
    expansion = j["expansion"].get<Expansion>();
    R = j["F/Mp"];
    meff2 = j["meff2"];


    N = j["N"];
    L = j["L"];
    num_fields  = j["num_fields"];
    num_threads = j["num_threads"];

    dt = j["dt"];
    output_step = j["output_step"];
    total_step  = j["total_step"];

    precision = j["precision"];


    rnd = j["rnd"];
    c = j["c"];
    t0 = j["t0"];
    ini_amp = j["initial_amplitude"].get<std::vector<double>>();
    ini_vel = j["initial_velocity"].get<std::vector<double>>();


    lambda = j["lambda"];


    dx = 1.* L/N;
    a   = 1;
    da  = 0;
    t   = t0;
}