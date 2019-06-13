#ifndef _ALLOCATOR_H_
#define _ALLOCATOR_H_


#include <cmath>
#include "parameter.hpp"



template <typename dataType>
dataType** allocateData( int num_fields, std::size_t N, int dimension )
{
    if( dimension < 1 or dimension > 3 )
    {
        std::cerr << "Error: Simulation dimension must be 1, 2, or 3" << std::endl;
        exit(1);
    }

    std::size_t grid_size = pow(N,dimension);

    dataType** data = new dataType* [num_fields];
    data[0] = new dataType [num_fields*grid_size]();
    for( int n = 0; n < num_fields; ++n ) data[n] = data[0] + n*grid_size;

    return data;
}


template <typename dataType>
void deleteData( dataType**& data )
{
    delete [] data[0];
    delete [] data;
}



#endif