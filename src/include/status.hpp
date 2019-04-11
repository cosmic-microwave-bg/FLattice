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
        std::shared_ptr<CalculateEnergy> _energy;
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
        Status( double** f, std::shared_ptr<CalculateEnergy> calculate )
        : _time(0), _ofs("../status.txt", std::ios::trunc), _f(f), _energy(calculate), 
        _field_average(num_fields), _field_variance(num_fields), _energy_average(num_fields)
        {
            if( !_ofs ){
                std::cerr << "Failed to open 'status.txt' file !" << std::endl;
                exit(1);
            }else{
                _ofs << "   t ";
                if( expansion != Expansion::no_expansion ) _ofs << "   a ";
                for( int n = 0; n < num_fields; ++n ) _ofs << "field_ave["  << n << "] ";
                for( int n = 0; n < num_fields; ++n ) _ofs << "field_var["  << n << "] ";
                for( int n = 0; n < num_fields; ++n ) _ofs << "energy_ave[" << n << "] ";
                _ofs << "total_energy_ave" << std::endl;
            }
        }

        void update( int now=0  )
        {
            _time = now;
            updateFieldAverage     ();
            updateFieldVariance    ();
            updateCalculateAverage ();
        }

        void write( int now=0 )
        {
            if( _time != now ) std::cerr << "Status are not updated when you write !" << std::endl;
            
            double a = LeapFrogExpansion::_a;

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
        }
};


inline void write_VTK( double* f, std::string str, int loop )
{
	unsigned int size;
	std::stringstream ss;
	std::ofstream fout;
	
	if( DIMENSION == 1 ){
		ss << "../data/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".txt";
    	fout.open( ss.str().c_str() );	
 
    	for( int j = 0; j < N; j++ ){
			int idx = j;
			fout << idx*dx << " " << f[idx] << std::endl;
		}
    }else{
    	ss << "../data/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".vti";
    	fout.open( ss.str().c_str() );
    
  	  	fout << "<?xml version=\"1.0\"?>" << std::endl;
  		fout << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl <<std::endl;
    	switch( DIMENSION ){
    		case 2:
				size = sizeof(double) * pow(N, 2);
				fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
				fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\">" << std::endl;
    			break;
			case 3:
				size = sizeof(double) * pow(N, 3);
				fout << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 <<" 0 " << N-1 << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
				fout << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << N-1 << "\">" << std::endl;
				break;
    	}
    	fout << "<PointData Scalars=\"energy\">" << std::endl;
    	fout << "<DataArray type=\"Float64\" Name=\"energy\" format=\"appended\" offset=\"0\" />" << std::endl;
    	fout << "</PointData>" << std::endl;
	    fout << "<CellData>" << std::endl;
	    fout << "</CellData>" << std::endl;
	    fout << "</Piece>" << std::endl << std::endl;
	    
	    fout << "</ImageData>" << std::endl << std::endl;
	    fout << "<AppendedData encoding=\"raw\">" << std::endl;
	    fout << "_" ;
	    fout.close();
	    
	    fout.open( ss.str().c_str(), std::ios::binary | std::ios::app);
	    fout.write( (char*) &size, sizeof(unsigned int) );
	    fout.write( (char*) f, size );
		
	    fout.close();	
	    fout.open( ss.str().c_str(), std::ios::app);
	    fout << std::endl << "</AppendedData>" << std::endl;
		fout << "</VTKFile>" ;
	}
	fout.close();
}



#endif