#ifndef _PHYSICAL_QUANTITY_H_
#define _PHYSICAL_QUANTITY_H_

#include "allocator.hpp"
#include "calculation.hpp"



class PhysicalQuantity
{
    protected:
        std::string _name;
        double** _data;
        double*  _data_tot;

    public:
        PhysicalQuantity           (std::string name, int num_fields, int N, int dimension)
        : _name(name)
        {
            _data = allocateData<double>( num_fields, N, dimension );
            if(dimension == 1) _data_tot = new double [N]();
            if(dimension == 2) _data_tot = new double [N*N]();
            if(dimension == 3) _data_tot = new double [N*N*N]();
        }

        virtual     ~PhysicalQuantity ()
        { 
            deleteData(_data);
            delete[] _data_tot; 
        }

        virtual void calculate     ( double**f, double** df ) = 0;

        std::string name   ()      { return _name; }
        double** data      ()      { return _data; }
        double*  data      (int n) { return _data[n]; }
        double*  data_tot  ()      { return _data_tot; }
        double   average   (int n) { return calculateAverage<double>(_data[n], N, DIMENSION); }
        double   average   ()      { return calculateAverage<double>(_data_tot, N, DIMENSION); }
        double   variance  (int n) { return calculateVariance<double>(_data[n], average(N), N, DIMENSION); }
        double   variance  ()      { return calculateVariance<double>(_data_tot, average(), N, DIMENSION); }
};



#endif