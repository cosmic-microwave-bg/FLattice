#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <sstream>
#include <fftw3.h>
#include "parameter.hpp"
#include "utilities.hpp"


void DFT_c2r( double** f )
{
	  fftw_plan p;
		double* out;
    fftw_complex *in;
		size_t in_size;
	  
		switch( dim )
		{
			case 1:
				in_size = sizeof(fftw_complex) * (N/2+1);
				in  = (fftw_complex*)fftw_malloc( in_size );
				out = new double [N]();
				p = fftw_plan_dft_c2r_1d( N, in, out, FFTW_ESTIMATE );
				break;
			case 2:
				in_size = sizeof(fftw_complex) * N * (N/2+1);
				in  = (fftw_complex*)fftw_malloc( in_size );
				out = new double [N*N]();
				p = fftw_plan_dft_c2r_2d( N, N, in, out, FFTW_ESTIMATE );
				break;
			case 3:
				in_size = sizeof(fftw_complex) * N * N * (N/2+1);
				in  = (fftw_complex*)fftw_malloc( in_size );
				out = new double [N*N*N]();
				p = fftw_plan_dft_c2r_3d( N, N, N, in, out, FFTW_ESTIMATE );
				break;
		}

		for( int i = 0; i < num_fields; ++i )
		{
				std::mt19937 mt( rnd );
				std::uniform_real_distribution<> rand( 0, 2*M_PI );
				
				// Create input data
				switch( dim ){
					case 1:
						#pragma omp parallel for schedule( static ) num_threads( num_threads )
						for( int j = 0; j < N/2+1; ++j ){
							int idx = j;
							double phase = rand(mt);
							in[idx][0] = f[i][idx] * cos( phase );
							in[idx][1] = f[i][idx] * sin( phase );
						}
						break;
					case 2:
						#pragma omp parallel for schedule( static ) num_threads( num_threads )
						for( int j = 0; j < N; ++j ){
							for( int k = 0; k < N/2+1; ++k ){
									int idx = j*N + k;
									double phase = rand(mt);
									in[idx][0] = f[i][idx] * cos( phase );
									in[idx][1] = f[i][idx] * sin( phase );
							}
						}
						break;
					case 3:
						#pragma omp parallel for schedule( static ) num_threads( num_threads )
						for( int j = 0; j < N; ++j ){
							for( int k = 0; k < N; ++k ){
									for( int l = 0; l < N/2+1; ++l ){
											int idx = (j*N + k)*N + l;
											double phase = rand(mt);
											in[idx][0] = f[i][idx] * cos( phase );
											in[idx][1] = f[i][idx] * sin( phase );
									}
							}
						}
						break;
				}
		
				fftw_execute(p);
	
				// Set output data
        #pragma omp parallel for schedule( static ) num_threads( num_threads )
        for( int j = 0; j < N; ++j ){
            switch( dim ){
							case 1:
                int idx = j;
                f[i][idx] = out[idx]/N;
								break;
            	case 2:
                for( int k = 0; k < N; ++k ){
                    int idx = j*N + k;
                    f[i][idx] = out[idx]/(N*N);
                }
								break;
            	case 3:
                for( int k = 0; k < N; ++k ){
                    for( int l = 0; l < N; ++l ){
                        int idx = (j*N + k)*N + l;
                        f[i][idx] = out[idx]/(N*N*N);
                    }
                }
								break;
						}
        }
    }

		if( p ) fftw_destroy_plan(p);
		if( in ) fftw_free(in);
		delete[] out;
}



void write_VTK( double* f, std::string str, int loop )
{
	unsigned int size;
	std::stringstream ss;
	std::ofstream fout;
	
	if( dim == 1 ){
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
		ofs.open( "../status.txt", std::ios::trunc );

		ofs << std::setw(3) << std::right << "  t ";
		if( expansion ) ofs << "  a ";
		for( int i = 0; i < num_fields; ++i ) ofs << "field_ave["  << i << "] ";
		for( int i = 0; i < num_fields; ++i ) ofs << "field_var["  << i << "] ";
		for( int i = 0; i < num_fields; ++i ) ofs << "energy_ave[" << i << "] ";
		ofs << "total_energy_ave" << std::endl;
	}
	else ofs.open( "../status.txt", std::ios::app );
	
	ofs << std::setw(3) << std::right << t << " ";
	if( expansion )
	{
		ofs << std::setw(3) << std::right << a << " ";
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


