project(vipss)
cmake_minimum_required(VERSION 2.8)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

# Shayan: this should already work out of the box...
SET(NLOPT_INCLUDE_DIRS "/usr/local/include/")
SET(NLOPT_LIB_DIR "/usr/local/")
SET(NLOPT_LIB nlopt)

# OLD
#SET(ARMADILLO_INCLUDE_DIRS "/Users/Research/Geometry/RBF/external/armadillo/include/")
#SET(ARMADILLO_LIB_DIRS "/Users/Research/Geometry/RBF/external/armadillo/")
#SET(ARMADILLO_LIB armadillo BLAS LAPACK)
# Shayan
find_package(Armadillo REQUIRED)
set(ARMADILLO_LIB_DIRS "")
set(ARMADILLO_LIB ${ARMADILLO_LIBRARIES})

include_directories(${NLOPT_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ./src/surfacer)
aux_source_directory(. MAIN)
aux_source_directory(./src SRC_LIST)
aux_source_directory(./src/surfacer SURFACER_LIST)

link_directories(${ARMADILLO_LIB_DIRS} ${NLOPT_LIB_DIR})
add_executable(${PROJECT_NAME} ${SRC_LIST} ${MAIN} ${SURFACER_LIST})

target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIB} ${NLOPT_LIB})
