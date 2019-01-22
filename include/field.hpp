#ifndef _FIELD_H_
#define _FIELD_H_

#include <cmath>
#include "parameter.hpp"


void initialize( double**& f, double**& df );
void finalize( double**& f, double**& df );

class Field
{
	private:
		double _dx;
		int _num_fields;
		int _num_threads;
        int _rnd;
		double* _average;
		double* _variance;

	public:
		Field ( double dx, int num_fields, int num_threads, int rnd );
        
		double gradient_energy  ( double* f ) const;
		double potential_energy ( double* f, double a ) const;
		double gradient  ( double* f, int d, int j, int k, int l ) const;
		double laplacian ( double* f, int j, int k, int l ) const;

		double V   ( double** f, int i, int idx ) const 
		{ return ( 1 + K * log( (pow(f[0][idx], 2) + pow(f[1][idx], 2))/(2*MG*MG) ) ) * (pow(f[0][idx], 2) + pow(f[1][idx], 2))/2; }
		double dV  ( double** f, int i, int idx ) const 
		{ return ( 1 + K * pow(f[i][idx], 2) / (pow(f[0][idx], 2) + pow(f[1][idx], 2)) + K * log( (pow(f[0][idx], 2) + pow(f[1][idx], 2))/(2*MG*MG) ) ) * f[i][idx]; }
		double aV  ( double** f, double a, int i, int idx ) const 
		{ return ( 1 + K * log( (pow(f[0][idx], 2) + pow(f[1][idx], 2))/(2*MG*MG*a*a) ) ) * pow(f[i][idx], 2)/(2*a*a); }
		double adV ( double** f, double a, int i, int idx )
		{ return ( 1 + K * pow(f[i][idx], 2) / (pow(f[0][idx], 2) + pow(f[1][idx], 2)) + K * log( (pow(f[0][idx], 2) + pow(f[1][idx], 2))/(2*a*a*MG*MG) ) ) * f[i][idx]/a; }
		
		/*double adV1 ( double** f, double a, int i, int idx )
		{ return ( 1 + K * pow(f[1][idx], 2) / (pow(f[0][idx], 2) + pow(f[1][idx], 2)) + K * log( (pow(f[0][idx], 2) + pow(f[1][idx], 2))/(2*a*a*MG*MG) ) ) * f[1][idx]/a; }

		double (Field::*adV[2]) ( double**, double, int, int ) = { &Field::adV0, &Field::adV1 };*/

		double average  ( double* f, int i );
		double variance ( double* f, int i );
};



#endif