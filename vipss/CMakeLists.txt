project(vipss)
cmake_minimum_required(VERSION 2.8)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 ")

find_package(nlopt REQUIRED)
set(NLOPT_LIB_DIR "")
set(NLOPT_LIB ${NLOPT_LIBRARIES})

find_package(Armadillo REQUIRED)
set(ARMADILLO_LIB_DIRS "")
set(ARMADILLO_LIB ${ARMADILLO_LIBRARIES})


include_directories(${NLOPT_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ./src/surfacer)
aux_source_directory(. MAIN)
aux_source_directory(./src SRC_LIST)
aux_source_directory(./src/surfacer SURFACER_LIST)

LINK_DIRECTORIES(${ARMADILLO_LIB_DIRS} ${NLOPT_LIB_DIR})
add_executable(${PROJECT_NAME} ${SRC_LIST} ${MAIN} ${SURFACER_LIST})

target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIB} ${NLOPT_LIB})
