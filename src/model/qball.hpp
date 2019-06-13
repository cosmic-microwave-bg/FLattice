#ifndef _QBALL_H_
#define _QBALL_H_

#include "model.hpp"


class Qball final: public Model
{
    public:
        Qball( std::string name ): Model(name) {}
        
        double V   ( double** f, int n, int idx, double a=1 ) const override
        { return ( 1 + K*log( (pow(f[0][idx],2) + pow(f[1][idx],2)) / (2*a*a)) ) * pow(f[n][idx],2)/(2*a*a); }
        double dV  ( double** f, int n, int idx, double a=1 ) const override
        { return ( 1+K + K*log( (pow(f[0][idx],2) + pow(f[1][idx],2)) / (2*a*a)) ) * f[n][idx]/a; }
        
};



#endif