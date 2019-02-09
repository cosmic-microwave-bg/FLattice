cmake_minimum_required ( VERSION 3.10 )

set ( CMAKE_CXX_COMPILER icpc )
set ( CMAKE_CXX_STANDARD 17 )
set ( CMAKE_CXX_FLAGS "-fast -qopenmp -Wall" )

project( FLattice CXX )

# Add include files (-I option)
include_directories ( ${PROJECT_SOURCE_DIR}/include )
include_directories ( /opt/fftw/include )

# Specify the Library directory (-L option)
link_directories ( /opt/fftw/lib )

# Add executing files
file ( GLOB lib_codes ${PROJECT_SOURCE_DIR}/lib/*.cpp )
add_executable ( ${PROJECT_NAME} main.cpp ${lib_codes} )

# Link external libraries
target_link_libraries ( ${PROJECT_NAME} fftw3 )


# Link project "library" when compile FLattice
# target_link_libraries ( FLattice library_code )

# Add sub-directory
# add_subdirectory ( lib )