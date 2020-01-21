file (GLOB SRC *.cpp)

add_library (base STATIC ${SRC})
set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fPIC")