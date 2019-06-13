#ifndef _LEAPFROG_IN_MINKOWSKY_H_
#define _LEAPFROG_IN_MINKOWSKY_H_

#include "leapfrog_base.hpp"



class LeapFrogInMinkowsky final: public LeapFrogBase
{
    public:
        LeapFrogInMinkowsky ( std::shared_ptr<Model> model, int precision )
        : LeapFrogBase(model, precision) {}
        
        void evolution  ( double** f, double** df, int output_step ) override;
};


class LeapFrogWithABC final: public LeapFrogBase
{
    private:
        double sphericalSymLaplacian  ( double* f, int i );
        void   evolFieldDerivsWithABC ( double **f, double **df, const double h );

    public:
        LeapFrogWithABC( std::shared_ptr<Model> model, int precision )
        : LeapFrogBase(model, precision)
        {
            if( DIMENSION != 1 ){
                std::cerr << "DIMENSION must be 1 when you use spherical symmetry." << std::endl;
                exit(1);
            }
        };

        void evolution   ( double** f, double** df, int output_step ) override;
};



#endif