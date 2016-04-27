project (FractureRB)
cmake_minimum_required (VERSION 2.8)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX ".." CACHE PATH "comment" FORCE)
endif()

option(BUILD_HYENA "If on, HyENA library parts are built from source, otherwise link to a pre-built version" OFF)
if( BUILD_HYENA )
	set(LIBTYPE SHARED)
endif()

add_definitions(-DASSEMBLE_GALERKIN_OMP) # use OpenMP in HyENA to assemble SGBEM matrices
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
	if( BUILD_HYENA )
		SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fopenmp")
	endif()
endif()

# my Header-only libraries are all in one directory, find that dir by searching for Eigen/Eigen
find_path(HLIB_INC Eigen/Eigen)
include_directories (${HLIB_INC}) # takes care of other header-only libs, need Eigen, TCLAP, and parts of Boost

# the HyENA library is in the project dir, because we've made a few adjustments to it
include_directories (${PROJECT_SOURCE_DIR})
# and here's where all the other headers are
include_directories (${PROJECT_SOURCE_DIR}/src)

# OpenVDB libraries and includes
find_path(OpenVDBinclude openvdb/openvdb.h)
include_directories(${OpenVDBinclude})
add_definitions(-DOPENVDB_STATICLIB) # use static OpenVDB lib
find_library(OpenVDBlib libOpenVDB)
find_library(Halflib Half)
find_library(tbblib tbb)
find_library(zlib zlib)
set(OpenVDB_libs ${OpenVDBlib} ${Halflib} ${tbblib} ${zlib})

# Bullet libraries and includes
add_definitions(-DBT_USE_DOUBLE_PRECISION)
find_path(BULLET_INCLUDE bullet/btBulletDynamicsCommon.h)
include_directories (${BULLET_INCLUDE}/bullet)
# Bullet core libraries
find_library(BULLET_DYNAMICS BulletDynamics)
find_library(BULLET_COLLISION BulletCollision)
find_library(BULLET_LINMATH LinearMath)
set ( BULLET_LIBS ${BULLET_DYNAMICS} ${BULLET_COLLISION} ${BULLET_LINMATH} ) # order is very important for GCC builds! (MSVC is ok either way)
# Bullet Extras libs
find_library(BULLET_IMPORTER BulletWorldImporter)
find_library(BULLET_LOADER BulletFileLoader)
set ( BULLET_EXT_LIBS  ${BULLET_IMPORTER} ${BULLET_LOADER} )


if( BUILD_HYENA )
  add_library ( HyENAlib2 ${LIBTYPE}
	"hyena/wrap/HyENAWrapper.cpp"
	"hyena/wrap/QuadratureWrapper.cpp"
  )
  install(TARGETS HyENAlib2 RUNTIME DESTINATION bin LIBRARY DESTINATION bin ARCHIVE DESTINATION bin) # throw this in the same dir as the binaries
else()
  find_library ( HYENA_LIB2 HyENAlib2 )
endif()

add_library ( FractureLib STATIC
	"src/FractureBEM.cpp"
	"src/FractureModel.cpp"
	"src/Reader.cpp"
	"src/Reader_VCG.cpp"
	"src/PostProcessor.cpp"
	"src/VDBWrapper.cpp"
	"src/VDBWrapper_mesh.cpp"
	"src/SubsampledCrackTip.cpp"
	"src/MaterialModel.cpp"
	"vcg/wrap/ply/plylib.cpp"
)
if( BUILD_HYENA )
  target_link_libraries (FractureLib HyENAlib2 )
else()
  target_link_libraries (FractureLib ${HYENA_LIB2} )
endif()
install(TARGETS FractureLib RUNTIME DESTINATION bin LIBRARY DESTINATION bin ARCHIVE DESTINATION bin) # throw this in the same dir as the binaries

FILE(GLOB SRC_HEADERS ${PROJECT_SOURCE_DIR}/src/*.h) # for convenience -> puts the headers into the MSVC project list
add_executable (FractureRB ${SRC_HEADERS}
	"src/main.cpp"
	"src/BulletWrapper.cpp"
	"src/FractureRB.cpp"
	"src/FractureRB_fragments.cpp"
	"src/ColliderData.cpp"
)
target_link_libraries (FractureRB  FractureLib ${OpenVDB_libs} ${BULLET_EXT_LIBS} ${BULLET_LIBS} ) # order is important for GCC
install(TARGETS FractureRB RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

add_executable (FractureRB_vis ${SRC_HEADERS}
	"src/mainVis.cpp"
	"src/VDBLoader.cpp"
)
target_link_libraries (FractureRB_vis FractureLib ${OpenVDB_libs} )
install(TARGETS FractureRB_vis RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

# old FractureBEM build targets -- mostly used for testing
#
add_executable (FractureBEM ${SRC_HEADERS}
	"src/mainFractureBEM.cpp"
)
target_link_libraries (FractureBEM FractureLib ${OpenVDB_libs} )
install(TARGETS FractureBEM RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
#
add_executable (SegmentVDB
	"src/mainSegment.cpp"
	"src/mySegment.h"
)
target_link_libraries (SegmentVDB ${OpenVDB_libs})
install(TARGETS SegmentVDB RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
#
add_executable (FractureBEM_vis ${SRC_HEADERS}
	"src/mainVisFractureBEM.cpp"
	"src/VDBLoader.cpp"
)
target_link_libraries (FractureBEM_vis FractureLib ${OpenVDB_libs})
install(TARGETS FractureBEM_vis RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
