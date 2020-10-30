/**
 * @file parameter.cpp
 * @brief See the explanation of the coresponding header file.
 */
#include <typeinfo>
#include "parameter.hpp"

#include <string.h>//file nameの比較とかに必要
#include <sstream>//string stream
#include <sys/stat.h>//mkdir

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

//global file name
std::string outfoldername;//出力フォルダの名前を決定。output_outputも取り扱えるように組む outfoldername.c_str()
std::string statusname;

//box parameters


double dx;
double a, da;
double t;


//model paramters 
double lambda;

void setParameters(char* filename)
{
    nlohmann::json j;
    std::stringstream inputname;
    std::stringstream outputfoldername;
    if(strcmp(filename, "")==0){
      inputname <<  "../src/parameter/parameter.json";
      outputfoldername<< "../output" << filename  ;
    }else{
  		inputname << "../" << filename  << ".json";
      outputfoldername<< "../output_" << filename  ;
    }
    mkdir(outputfoldername.str().c_str(), 0755);//filenameに対応するファイル名の設定
    outfoldername = outputfoldername.str();
    std::ifstream input( inputname.str().c_str() );
    if (input.fail()){ std::cout<<"We can not find the json"<<std::endl;  exit(0);  }
    input >> j;

    //面倒だが、いい定義が思いつかない。
    std::stringstream ssStatusname;
    ssStatusname << outfoldername << "/status.txt";
    statusname = ssStatusname.str();


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