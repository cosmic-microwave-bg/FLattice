#ifndef _ENERGY_H_
#define _ENERGY_H_

#include "physical_quantity.hpp"


class Energy final: public PhysicalQuantity
{
    public:
        Energy (std::string name, int num_fields, int N, int dimension)
        : PhysicalQuantity(name, num_fields, N, dimension) {}

        void calculate (double** f, double** df) override;
};


#endif