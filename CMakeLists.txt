cmake_minimum_required(VERSION 3.9)
project(RADADMM)

set(CMAKE_CXX_STANDARD 11)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
#INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include /data/home/xjsjleiyongmei/xjy/local/include)
#LINK_DIRECTORIES(/data/home/xjsjleiyongmei/xjy/local/lib)
link_libraries(pthread zmq)

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_FILE)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/optimizer SRC_FILE)
aux_source_directory(${PROJECT_SOURCE_DIR}/src/utils SRC_FILE)

add_executable(radadmm ${SRC_FILE})
