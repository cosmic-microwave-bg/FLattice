/**
 * @file calculation.hpp
 * @brief The classes to calculate pysical values, such as energy, charge, and so on, are implemented.
 */
#ifndef _CALCULATION_H_
#define _CALCULATION_H_

#include <memory>
#include "model.hpp"
#include "evolution.hpp"



class CalculateBase
{
    protected:
        double** _data;
        double*  _total;
        double*  _average;

    public:
        CalculateBase(): _average(new double [num_fields]())
        {
            _data = allocateData<double>( num_fields, N, DIMENSION );
            if( DIMENSION == 1 ) _total = new double [N]();
            if( DIMENSION == 2 ) _total = new double [N*N]();
            if( DIMENSION == 3 ) _total = new double [N*N*N]();
        }
};


class CalculateEnergy: public CalculateBase
{
    private:
    public:
        CalculateEnergy( std::shared_ptr<Model> model, double** f, double** df );
        double gradient_energy( double* f );
        //writeVTI();
};



#endif