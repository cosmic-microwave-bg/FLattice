cmake_minimum_required ( VERSION 3.10 )

set ( CMAKE_CXX_COMPILER icpc )
set ( CMAKE_CXX_STANDARD 17 )
set ( CMAKE_CXX_FLAGS "-ipo -xhost -no-prec-div -qopenmp -parallel" )  # -guide


project( FLattice CXX )


# Add include files (-I option)
# Include files must be added before call "add_executable"
include_directories ( ${PROJECT_SOURCE_DIR}/include )
include_directories ( ${PROJECT_SOURCE_DIR}/src/include )
include_directories ( ${PROJECT_SOURCE_DIR}/src/evolution_shceme/include )
include_directories ( ${PROJECT_SOURCE_DIR}/src/physical_quantity/include )
include_directories ( ${PROJECT_SOURCE_DIR}/src/model )
include_directories ( ${PROJECT_SOURCE_DIR}/src/parameter )

# Specify the Library directory (-L option)
link_directories ( /usr/local/lib )

# Add executing files
file ( GLOB src_codes ${PROJECT_SOURCE_DIR}/src/*.cpp )
file ( GLOB evol_codes ${PROJECT_SOURCE_DIR}/src/evolution_shceme/*.cpp )
file ( GLOB phys_codes ${PROJECT_SOURCE_DIR}/src/physical_quantity/*.cpp )
file ( GLOB para_codes ${PROJECT_SOURCE_DIR}/src/parameter/*.cpp )
add_executable ( ${PROJECT_NAME} ${src_codes} ${evol_codes} ${phys_codes} ${para_codes} )

# Link external libraries
# target_link_libraries ( ${PROJECT_NAME} fftw3 vtkCommonDataModel-8.2.1 vtkCommonCore-8.2.1 vtkIOXML-8.2.1 )
target_link_libraries ( ${PROJECT_NAME} fftw3 )


# Add sub-directory
# add_subdirectory ( lib )