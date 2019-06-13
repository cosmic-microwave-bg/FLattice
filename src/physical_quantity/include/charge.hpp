#ifndef _CHARGE_H_
#define _CHARGE_H_

#include <cmath>
#include <queue>
#include <vector>
#include <utility>
#include <fstream>
#include "physical_quantity.hpp"



class Charge: public PhysicalQuantity
{
    private:
        double _Qini;

    public:
        Charge( std::string name, int num_fields, int N, int dimension, double**f, double** df )
        : PhysicalQuantity(name, num_fields, N, dimension)
        {
            calculate(f, df);
            _Qini = calculateAverage(_data_tot, N, DIMENSION);
        }

        void calculate( double**f, double** df) override;
};



#endif