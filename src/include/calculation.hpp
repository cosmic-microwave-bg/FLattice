/**
 * @file calculation.hpp
 * @brief The classes to calculate pysical values, such as energy, charge, and so on, are implemented.
 */
#ifndef _CALCULATION_H_
#define _CALCULATION_H_

#include <memory>
#include "model.hpp"
#include "evolution.hpp"
#include "write.hpp"



template <typename dataType>
dataType calculateAverage( dataType* data, int N, int dimension )
{
    dataType average = 0;

    #pragma omp parallel for reduction(+:average) schedule( static ) num_threads ( num_threads )
    for( int i = 0; i < N; ++i ){
		switch( dimension )
		{
			case 1:
				int idx = i;
				average += data[idx];
				break;
			case 2:
				for( int j = 0; j < N; ++j ){
					int idx = i*N+j;
					average += data[idx];
				}
				break;
			case 3:
				for( int j = 0; j < N; ++j )
					for( int k = 0; k < N; ++k ){
						int idx = (i*N+j)*N+k;
						average += data[idx];
					}
				break;
		}
	}
	for( int i = 0; i < dimension; ++i ) average /= N;
	
	return average;
}


template <typename dataType>
dataType calculateVariance( dataType* data, dataType average, int N, int dimension )
{
    dataType variance = 0;

    #pragma omp parallel for reduction(+:variance) schedule( static ) num_threads ( num_threads )
    for( int i = 0; i < N; ++i ){
		switch( dimension )
        {
			case 1:
				int idx = i;
				variance += pow( data[idx]-average, 2 );
				break;
			case 2:
				for( int j = 0; j < N; ++j ){
					int idx = i*N+j;
    				variance += pow( data[idx]-average, 2 );
				}
				break;
			case 3:
				for( int j = 0; j < N; ++j ){
					for( int k = 0; k < N; ++k ){
						int idx = (i*N+j)*N+k;
    				    variance += pow( data[idx]-average, 2 );
					}
				}
				break;
		}
	}
    for( int i = 0; i < dimension; ++i ) variance /= N;
	
	return variance;
}


/** 
 * @class CalculateBase
 * @brief Abstract base class to calculate phycal quantites.
 */
class CalculateBase
{
    protected:
		double** _f;
		double** _df;
        double** _data;      //!< The physical quantity you want to calculate is stored.
        double*  _data_tot;  //!< Total of _data at each grid point.

		void         calculate     ();
		/**
		 * @brief Pure virtual function to be implemented in the sub class.
		 */
		virtual void calculateData ( int n, int i, int j=0, int k=0 ) = 0;

    public:
        CalculateBase              ( double** f, double** df ): _f(f), _df(df)
        {
			_data = allocateData<double>( num_fields, N, DIMENSION );
            if( DIMENSION == 1 ) _data_tot = new double [N]();
            if( DIMENSION == 2 ) _data_tot = new double [N*N]();
            if( DIMENSION == 3 ) _data_tot = new double [N*N*N]();
        }
        virtual     ~CalculateBase ()
		{ 
			deleteData(_data);
			delete[] _data_tot; 
		}

		void        update         ()        { calculate(); }
		double      average        ( int n ) { return calculateAverage<double>(_data[n], N, DIMENSION); }
};


class CalculateEnergy: public CalculateBase
{
    private:
		std::shared_ptr<Model> _model;
		double      gradientEnergy ( const double* f, int i, int j=0, int k=0 ) const;

    public:
        CalculateEnergy            ( std::shared_ptr<Model> model, double** f, double** df )
		:CalculateBase(f, df), _model(model) { }
        
		void        calculateData  ( int n, int i, int j=0, int k=0 ) override;
		void        write          ( int loop )
		{
			writeVTI<double>( _data, num_fields, "energy", loop);
			writeVTI<double>( _data_tot, "total_energy", loop  );
		}
};



#endif