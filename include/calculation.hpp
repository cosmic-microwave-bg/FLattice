#ifndef _CALCULATION_H_
#define _CALCULATION_H_

//#include <valarray>

#include "field.hpp"
#include "evolution.hpp"


class Calculation
{
    protected:
        int _num_fields;
        double* _average;
        double* _variance;

    public:
        double** value;
        Calculation( int num_fields );

        double average( int i ){ return _average[i]; }
        double variance( int i ){ return _variance[i]; }
};


class Energy: public Calculation
{
    private:
        double _dx;
        
    public:
        Energy( Field* field, LeapFrog* leapfrog, double** f, double** df, int num_fields, double dx );
        double gradient_energy( double* f );
};


#endif