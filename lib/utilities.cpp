#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "parameter.hpp"
#include "utilities.hpp"


void write_VTK( double* f, std::string str, int loop )
{
	unsigned int size;
	std::stringstream ss;
	std::ofstream fout;
	
	if( dim == 1 ){
		ss << "./data/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".txt";
    	fout.open( ss.str().c_str() );	
 
    	for( int j = 0; j < N; j++ ){
			int idx = j;
			fout << idx*dx << " " << f[idx] << std::endl;
		}
    }else{
    	ss << "./data/" << str << "." << std::setw(4) << std::setfill('0') << loop+1 <<".vti";
    	fout.open( ss.str().c_str() );
    
  	  	fout << "<?xml version=\"1.0\"?>" << std::endl;
  		fout << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl <<std::endl;
    	switch( dim ){
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


void write_status( Field* field, LeapFrog* leapfrog, Energy* energy, double** f, double t )
{
	double a = leapfrog->a();
	std::ofstream ofs;
	
	if( t == t0 )
	{
		ofs.open( "status.txt", std::ios::trunc );
		ofs << "t ";
		if( expansion ) ofs << "a ";
		for( int i = 0; i < num_fields; ++i ) ofs << "field_average[" << i << "] ";
		for( int i = 0; i < num_fields; ++i ) ofs << "field_variance[" << i << "] ";
		for( int i = 0; i < num_fields; ++i ) ofs << "energy_average[" << i << "] ";
		ofs << "total_energy_average" << std::endl;
	}
	else ofs.open( "status.txt", std::ios::app );
	
	ofs << std::setprecision(4) << t << " ";
	if( expansion )
	{
		ofs << a << " ";
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->average(f[i], i)/a << " ";
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->variance(f[i], i)/(a*a) << " ";
	}
	else
	{
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->average(f[i], i) << " ";
		for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << field->variance(f[i], i) << " ";
	}
	for( int i = 0; i < num_fields; ++i ) ofs << std::showpos << std::scientific << std::setprecision(4) << energy->average(i) << " ";
	ofs << energy->total_average() << std::endl;
}


