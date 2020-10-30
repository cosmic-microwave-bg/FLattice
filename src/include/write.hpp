#ifndef _WRITE_H_
#define _WRITE_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include "parameter.hpp"


/** You can use VTK Library instead of the below code.
 *  However, VTK Library code fairly slows down the output process.
 */ 

template <typename dataType>
void writeVTI(dataType** f, int num_fields, std::string name, int loop)
{
	unsigned int size;
	std::stringstream filename;
	std::ofstream ofs;

    filename << outfoldername.c_str() << "/" << name+'_' << std::setw(4) << std::setfill('0') << loop;

	if( DIMENSION == 1 ){
		filename << ".txt";
    	ofs.open( filename.str().c_str() );	
 
    	for( int i = 0; i < N; ++i ){
            ofs << i*dx << " ";
            for( int n = 0; n < num_fields; ++n ) ofs << f[n][i] << " ";
            ofs << std::endl;
		}
    }else{
    	filename << ".vti";
    	ofs.open( filename.str().c_str() );
    
  	  	ofs << "<?xml version=\"1.0\"?>" << std::endl;
  		ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl; 
    	switch( DIMENSION )
        {
    		case 2:
                size = sizeof(dataType)*pow( N, 2 );
                ofs << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
                ofs << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\">" << std::endl;
                break;
            case 3:
                size = sizeof(dataType)*pow( N, 3 );
                ofs << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 <<" 0 " << N-1 << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
                ofs << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << N-1 << "\">" << std::endl;
                break;
    	}
    	ofs << "<PointData>" << std::endl;

        for( int n = 0; n < num_fields; ++n ){
            std::size_t offset = n * (size + sizeof(unsigned int));
            ofs << "<DataArray type=\"Float64\" Name=\"" << name << "[" << n << "]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\"/>" << std::endl;
        }
        ofs << "</PointData>" << std::endl;
        ofs << "<CellData>" << std::endl;
        ofs << "</CellData>" << std::endl;
        ofs << "</Piece>" << std::endl;
        ofs << "</ImageData>" << std::endl;

        ofs << "<AppendedData encoding=\"raw\">" << std::endl;
        ofs << "_" ;
        ofs.close();
        ofs.open( filename.str().c_str(), std::ios::binary | std::ios::app);
        for( int n = 0; n < num_fields; ++n ){
            ofs.write( (char*) &size, sizeof(unsigned int) );
            ofs.write( (char*) f[n], size );
        }
        ofs.close();

        ofs.open( filename.str().c_str(), std::ios::app);
        ofs << "</VTKFile>" ;
    }
}


template <typename dataType>
void writeVTI(dataType* f, std::string name, int loop)
{
	unsigned int size;
	std::stringstream filename;
	std::ofstream ofs;

     filename << outfoldername.c_str() << "/" << name+'_' << std::setw(4) << std::setfill('0') << loop;

	if( DIMENSION == 1 ){
		filename << ".txt";
    	ofs.open( filename.str().c_str() );	
    	for( int i = 0; i < N; ++i ) ofs << i*dx << " " << f[i] << std::endl;
    }else{
    	filename << ".vti";
    	ofs.open( filename.str().c_str() );
    
  	  	ofs << "<?xml version=\"1.0\"?>" << std::endl;
  		ofs << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << std::endl; 
    	switch( DIMENSION )
        {
    		case 2:
                size = sizeof(dataType)*pow( N, 2 );
                ofs << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
                ofs << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << " 0 " << "\">" << std::endl;
                break;
            case 3:
                size = sizeof(dataType)*pow( N, 3 );
                ofs << "<ImageData WholeExtent=\"0 " << N-1 << " 0 " << N-1 <<" 0 " << N-1 << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dx << " " << dx << "\">" << std::endl;
                ofs << "<Piece Extent=\"0 " << N-1 << " 0 " << N-1 << " 0 " << N-1 << "\">" << std::endl;
                break;
    	}
    	ofs << "<PointData>" << std::endl;

        ofs << "<DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"1\" format=\"appended\" offset=\"0\"/>" << std::endl;
        ofs << "</PointData>" << std::endl;
        ofs << "<CellData>" << std::endl;
        ofs << "</CellData>" << std::endl;
        ofs << "</Piece>" << std::endl;
        ofs << "</ImageData>" << std::endl;

        ofs << "<AppendedData encoding=\"raw\">" << std::endl;
        ofs << "_" ;
        ofs.close();
        ofs.open( filename.str().c_str(), std::ios::binary | std::ios::app);
        ofs.write( (char*) &size, sizeof(unsigned int) );
        ofs.write( (char*) f, size );
        ofs.close();

        ofs.open( filename.str().c_str(), std::ios::app);
        ofs << "</VTKFile>" ;
    }
}



#endif