#ifndef _FIELD_H_
#define _FIELD_H_

#include "allocator.hpp"


struct Field
{
    double **f, **df;
    Field(int num_fields, int N, int dimension){
        f  = allocateData<double>(num_fields, N, DIMENSION);
        df = allocateData<double>(num_fields, N, DIMENSION);
    }

    ~Field() {
        deleteData(f);
        deleteData(df);
    }
};



#endif