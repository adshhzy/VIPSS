project(vipss)
cmake_minimum_required(VERSION 2.8)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 ")

SET(NLOPT_INCLUDE_DIRS "/usr/local/include/")
SET(NLOPT_LIB_DIR "/usr/local/")
SET(NLOPT_LIB nlopt)


#SET(SUITESPARSE_LIB_DIR "/usr/local/Cellar/suite-sparse/5.2.0/lib/")
#SET(SUITESPARSE_LIB cholmod colamd)

#SET(SUPERLU_LIB_DIR "/Users/Research/Geometry/RBF/external/superlu/")
#SET(SUPERLU_LIB superlu)

SET(ARMADILLO_INCLUDE_DIRS "/Users/Research/Geometry/RBF/external/armadillo/include/")
SET(ARMADILLO_LIB_DIRS "/Users/Research/Geometry/RBF/external/armadillo/")
SET(ARMADILLO_LIB armadillo BLAS LAPACK)

include_directories(${NLOPT_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ./src/surfacer)
aux_source_directory(. MAIN)
aux_source_directory(./src SRC_LIST)
aux_source_directory(./src/surfacer SURFACER_LIST)

#LINK_DIRECTORIES(${ARMADILLO_LIB_DIRS} ${SUITESPARSE_LIB_DIR} ${NLOPT_LIB_DIR} ${SUPERLU_LIB_DIR})
LINK_DIRECTORIES(${ARMADILLO_LIB_DIRS} ${NLOPT_LIB_DIR})
add_executable(${PROJECT_NAME} ${SRC_LIST} ${MAIN} ${SURFACER_LIST})

#target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIB} ${SUITESPARSE_LIB} ${NLOPT_LIB} ${SUPERLU_LIB})
target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIB} ${NLOPT_LIB})
