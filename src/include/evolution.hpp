/**
 * @file evolution.hpp
 * @brief The header file of classes necessary for the time evolution.
 */
#ifndef _EVOULUTION_H_
#define _EVOULUTION_H_

#include <memory>
#include "model.hpp"



/** 
 * @class Evolution
 * @brief The base class of time evolution of the fields.
 * 
 * We will implement the factory method to generate the evolution class ( LeapFrog, LeapFrogExpansion, ... )
 * by simply set the parameter 'expansion' and so on.
 */
class Evolution
{
	protected:
		std::shared_ptr<Model> _model;
		/** @brief           The function to calculate the laplacian by fourth-order centeral differencial scheme.
		 *  @param[in]  *f   The argument of the fields must be the pointer type for the vector optimization of the intel compiler!
		 */
		double   laplacian    ( const double* f, int i, int j=0, int k=0 ) const;
		void     evol_fields  ( double** f, double** df, const double h );  // This member function should not be here !

	public:
		explicit Evolution    ( std::shared_ptr<Model> model ): _model(model) {}
};


/** 
 * @class LeapFrog
 * @brief Evolve the fields by leapfrog method.
 */
class LeapFrog: public Evolution
{
	private:
		/** @brief          Evolve the conjugates of the fields df.
		 *  @param[in]  h   A constant depending on the time evolution precision.
		 */
        void     evol_field_derivs  ( double** f, double** df, const double h );
		
	public:
		LeapFrog                    ( std::shared_ptr<Model> model ): Evolution(model) 
		{ if( precision != 2 and precision != 4 ) std::cerr << "LeapFrog precision must be 2 or 4." << std::endl, exit(1); }
        void     evolution          ( double** f, double** df );
};


/** 
 * @class LeapFrogExpansion
 * @brief Evolve the fields by leapfrog method in the expanding universe.
 */
class LeapFrogExpansion: public Evolution
{
	private:
		void     evol_field_derivs   ( double** f, double** df, double h );
		void     evol_scale          ( double** f, double h );
		void     evol_scale_derivs   ( double h ) { _da += _dda * h*dt; }

	public:
		LeapFrogExpansion            ( std::shared_ptr<Model> model );
		void     evolution           ( double** f, double** df, double t );
		static double _a, _da, _dda;
};



#endif