cmake_minimum_required ( VERSION 3.10 )

set ( CMAKE_CXX_COMPILER icpc )
set ( CMAKE_CXX_STANDARD 17 )
set ( CMAKE_CXX_FLAGS "-fast -qopenmp -Wall" )
set ( INC_DIR ${PROJECT_SOURCE_DIR}/include/ )

project( FLattice CXX )

# Add executing files
file ( GLOB lib_codes ${PROJECT_SOURCE_DIR}/lib/*.cpp )
add_executable ( ${PROJECT_NAME} main.cpp ${lib_codes} )

include_directories ( ${PROJECT_SOURCE_DIR}/include )


# Link project "library" when compile FLattice
#target_link_libraries ( FLattice library_code )
# Add sub-directory
#add_subdirectory ( lib )