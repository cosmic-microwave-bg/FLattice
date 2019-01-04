#ifndef _FIELD_H_
#define _FIELD_H_


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
		Field( double dx, int num_fields, int num_threads, int rnd );
        
		double laplacian( double* f, int j, int k, int l ) const;
        double gradient( double* f, int d, int j, int k, int l ) const;
		double gradient_energy( double* f ) const;
		double potential_energy( double* f, double a ) const;

		double V(double** f, int i, int idx ) const { return f[i][idx]*f[i][idx]/2; } 
		double dV( double** f, int i, int idx ) const { return f[i][idx]; }
		double aV(double** f, double a, int i, int idx ) const { return f[i][idx]*f[i][idx]/(2*a*a); } 
		double adV( double** f, double a, int i, int idx ) const { return f[i][idx]/a; };

		double average( double* f, int i );
		double variance( double* f, int i );
};



#endif