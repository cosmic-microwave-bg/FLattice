/**
 * @file write.hpp
 * @brief The header only file of the data writer using the vtk library.
 * 
 * This file is not used for now because of the speed.
 */
#ifndef _WRITEVTK_H_
#define _WRITEVTK_H_

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <vtk-8.1/vtkSmartPointer.h>
#include <vtk-8.1/vtkImageData.h>
#include <vtk-8.1/vtkPointData.h>
#include <vtk-8.1/vtkXMLImageDataWriter.h>
#include <vtk-8.1/vtkDoubleArray.h>             // including these two headers prevent the intel compiler's optimization!
#include <vtk-8.1/vtkAOSDataArrayTemplate.h>
#include "parameter.hpp"



/**
 * @brief  Create the .vti file of 'data'.
 *
 * Including this file, that is, vtk library, may output '-ipo warning' when compile
 * by the intel compiler because of its bugs. 
 */
template <typename dataType>
void writeVTI ( dataType** const data, int num_fields, std::string name, int loop )
{
	std::vector<vtkSmartPointer<vtkAOSDataArrayTemplate<dataType>>> dataArray(num_fields);
	for( int n = 0; n < num_fields; ++n ){
		std::stringstream ss;
		ss << name << '[' << n << ']';
        
        dataArray[n] = vtkSmartPointer<vtkDoubleArray>::New();
		dataArray[n]->SetName( ss.str().c_str() );
		dataArray[n]->SetNumberOfComponents(1);
        switch( DIMENSION )
        {
            case 2:
                dataArray[n]->SetNumberOfTuples(N*N);
                dataArray[n]->SetArray(data[n], N*N, 1);
                break;
            case 3:
                dataArray[n]->SetNumberOfTuples(N*N*N);
                dataArray[n]->SetVoidArray(data[n], N*N*N, 1);
                break;
            default:
                std::cerr << "DIMESION must be 2 or 3 when you use function 'writeVTI'." << std::endl;
                exit(1);
        }
	}

    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    if( DIMENSION == 2 ){
        imageData->SetDimensions(N, N, 0);
        imageData->SetSpacing(dx, dx, 0);
    }
    if( DIMENSION == 3 ){
        imageData->SetDimensions(N, N, N);
        imageData->SetSpacing(dx, dx, dx);
    }
	for( int n = 0; n < num_fields; ++n ) imageData->GetPointData()->AddArray( dataArray[n] );

    std::stringstream ss;
	ss << "../data/" << name << std::setw(4) << std::setfill('0') << loop << ".vti";
	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName( ss.str().c_str() );
	writer->SetInputData( imageData );
	writer->SetDataModeToBinary();
	writer->Write();
}



#endif