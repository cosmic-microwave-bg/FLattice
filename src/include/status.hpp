/**
 * @file  status.hpp
 * @brief 
 */
#ifndef _STATUS_H_
#define _STATUS_H_

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <memory>
#include <valarray>
#include "model.hpp"
#include "calculation.hpp"
#include "write.hpp"


/**
 * @class Status
 * @brief Manage the field average, energy average and so on.
 * 
 * This class is temporary to observe the simulation correctness and may have bugs,
 * so use with caution.
 */
class Status
{
    private:
        int _time;  //!< To check wether the status is updeated when you call the member function 'write()'.
        std::ofstream _ofs;  //!< Instance to write the status in the status file.
        double** _f;   //!< Pointer to the fields f.
        std::unique_ptr<CalculateEnergy> _energy;
        std::valarray<double> _field_average;
        std::valarray<double> _field_variance;
        std::valarray<double> _energy_average;

        void updateFieldAverage()
        { for( int n = 0; n < num_fields; ++n ) _field_average[n] = calculateAverage<double>(_f[n], N, DIMENSION); }
        void updateFieldVariance()
        { for( int n = 0; n < num_fields; ++n ) _field_variance[n] = calculateVariance<double>(_f[n], _field_average[n], N, DIMENSION); }
        void updateCalculateAverage()
        { for( int n = 0; n < num_fields; ++n ) _energy_average[n] = _energy->average(n); }

    public:
        Status( double** f, std::unique_ptr<CalculateEnergy> energy )
        : _time(0), _ofs("../status.txt", std::ios::trunc), _f(f), _energy(std::move(energy)), 
        _field_average(num_fields), _field_variance(num_fields), _energy_average(num_fields)
        {
            if( !_ofs ){
                std::cerr << "Failed to open 'status.txt' file !" << std::endl;
                exit(1);
            }else{
                _ofs << "   t ";
                if( expansion != Expansion::no_expansion ) _ofs << "    a ";
                for( int n = 0; n < num_fields; ++n ) _ofs << "field_ave["  << n << "] ";
                for( int n = 0; n < num_fields; ++n ) _ofs << "field_var["  << n << "] ";
                for( int n = 0; n < num_fields; ++n ) _ofs << "energy_ave[" << n << "] ";
                _ofs << "total_energy_ave" << std::endl;
            }
        }

        void update( int now=0  )
        {
            _time = now;
            _energy->update();
            updateFieldAverage     ();
            updateFieldVariance    ();
            updateCalculateAverage ();
        }

        void write( int now=0 )
        {
            if( _time != now ) std::cerr << "Status are not updated when you write !" << std::endl;
            
            double a = LeapFrogExpansion::a;

            _ofs << std::noshowpos << std::fixed <<std::setprecision(2);
            _ofs << now*output_step*dt << "\t";
            
            if( expansion != Expansion::no_expansion )
            {
                _ofs << a << " ";
                _ofs << std::showpos << std::scientific << std::setprecision(3);
                for( int n = 0; n < num_fields; ++n ) _ofs << _field_average[n]/a      << "\t";
                for( int n = 0; n < num_fields; ++n ) _ofs << _field_variance[n]/(a*a) << "\t";
            }
            else
            {
                _ofs << std::showpos << std::scientific << std::setprecision(3);
                for( int n = 0; n < num_fields; ++n ) _ofs << _field_average[n]  << "\t";
                for( int n = 0; n < num_fields; ++n ) _ofs << _field_variance[n] << "\t";
            }
            for( int i = 0; i < num_fields; ++i ) _ofs << _energy_average[i] << "\t";
            _ofs << _energy_average.sum()/num_fields << std::endl;

            _energy->write(now);
            //writeVTI( _energy->data_sum, num_fields, "total_energy", loop);
        }
};



#endif