project (HyENAlib2)
cmake_minimum_required (VERSION 2.8)

#
# Use this CMakeLists file for building ONLY the HyENA lib
# -- for building both the lib and all of the FractureBEM binaries
# -- use the CMakeLists in the parent directory
#

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX ".." CACHE PATH "comment" FORCE)
endif()
set(LIBTYPE SHARED) # GCC (at least MinGW) won't properly link OpenMP into static libs

add_definitions(-DASSEMBLE_GALERKIN_OMP) # use OpenMP in HyENA to assembly SGBEM matrices
if(MSVC)
	# stuff for Visual Studio compiler
	add_definitions(-TP) # force MSVC to use C++ compiler, otherwise it falls back to ANSI C
	add_definitions(-D_USE_MATH_DEFINES) # make sure we have M_PI available
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp") # enable OpenMP
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0") # get rid of warnings
else()
	# stuff for GNU compiler
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp") # enable OpenMP
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w") # get rid of warnings
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fopenmp")
endif()

find_path(EIGEN_INC Eigen/Eigen)
include_directories (${EIGEN_INC}) # we also need some headers from boost from the same directory

# the HyENA library is in the main project dir
include_directories (${PROJECT_SOURCE_DIR}/..)
include_directories (${PROJECT_SOURCE_DIR}/../src)
add_library ( HyENAlib2 ${LIBTYPE}
 "wrap/HyENAWrapper.cpp"
 "wrap/QuadratureWrapper.cpp"
)

install(TARGETS HyENAlib2 RUNTIME DESTINATION bin LIBRARY DESTINATION bin ARCHIVE DESTINATION bin) # throw this in the same dir as the binaries
