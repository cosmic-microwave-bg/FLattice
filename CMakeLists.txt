cmake_minimum_required ( VERSION 3.10 )

set ( CMAKE_CXX_COMPILER icpc )
set ( CMAKE_CXX_STANDARD 17 )
#set ( CMAKE_CXX_FLAGS "-ipo -xhost -no-prec-div -qopenmp -Wall" )
set ( CMAKE_CXX_FLAGS "-ipo -xhost -no-prec-div -qopenmp" )


project( FLattice CXX )


# Add include files (-I option)
# Include files must be added before call "add_executable"
include_directories ( ${PROJECT_SOURCE_DIR}/src/include )

# Specify the Library directory (-L option)
link_directories ( /usr/local/lib )

# Add executing files
file ( GLOB src_codes ${PROJECT_SOURCE_DIR}/src/*.cpp )
add_executable ( ${PROJECT_NAME} ${src_codes} )

# Link external libraries
# target_link_libraries ( ${PROJECT_NAME} fftw3 vtkCommonDataModel-8.1.1 vtkCommonCore-8.1.1 )
target_link_libraries ( ${PROJECT_NAME} fftw3 )


# Add sub-directory
# add_subdirectory ( lib )