#ifndef _CREATE_PHYSICAL_QUANTITY_H_
#define _CREATE_PHYSICAL_QUANTITY_H_

#include <memory>
#include "parameter.hpp"
#include "energy.hpp"
#include "charge.hpp"


template <typename T=PhysicalQuantity>
std::unique_ptr<T> createPhysicalQuantity( std::string name, std::shared_ptr<Model> model, double** f, double** df )
{
    if( name == "energy" )
        return std::make_unique<Energy>(name, num_fields, N, DIMENSION, model);
    else if( name == "charge" )
        return std::make_unique<Charge>(name, num_fields, N, DIMENSION, f, df);
    
    return nullptr;
}



#endif