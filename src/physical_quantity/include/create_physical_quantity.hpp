#ifndef _CREATE_PHYSICAL_QUANTITY_H_
#define _CREATE_PHYSICAL_QUANTITY_H_

#include <memory>
#include "parameter.hpp"
#include "energy.hpp"
#include "charge.hpp"


template <typename T=PhysicalQuantity>
std::unique_ptr<T> createPhysicalQuantity(std::string name, double** f, double** df)
{
    if(name == "energy")
        return std::make_unique<Energy>(name, num_fields, N, DIMENSION);
    if(name == "charge"){
        if(num_fields != 2){
            std::cerr << "'num_fields' must be 2 when you set 'charge'." << std::endl;
            exit(1);
        }
        return std::make_unique<Charge>(name, num_fields, N, DIMENSION, f, df);
    }
    
    return nullptr;
}



#endif