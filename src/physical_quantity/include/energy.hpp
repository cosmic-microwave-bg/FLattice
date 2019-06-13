#ifndef _ENERGY_H_
#define _ENERGY_H_


#include "physical_quantity.hpp"
#include "model.hpp"


class Energy final: public PhysicalQuantity
{
    private:
        std::shared_ptr<Model> _model;

    public:
        Energy ( std::string name, int num_fields, int N, int dimension, std::shared_ptr<Model> model )
        : PhysicalQuantity(name, num_fields, N, dimension), _model(model) {}

        void calculate ( double** f, double** df ) override;
};



#endif