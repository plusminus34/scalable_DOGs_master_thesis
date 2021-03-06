cmake_minimum_required(VERSION 3.1)
project(dog_editor)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_CGAL            "Use CGAL"  ON)
option(LIBIGL_WITH_TRIANGLE            "Use Triangle"  ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)

set(CMAKE_BUILD_TYPE "Release")

find_package(LIBIGL REQUIRED QUIET)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -ffast-math -O2")

#set(EIGEN3_VERSION_OK TRUE)
#set(EIGEN3_INCLUDE_DIR "/Users/michaelrabinovich/libigl/external/eigen")
#find_package(ifopt REQUIRED)

#list(APPEND CMAKE_MODULE_PATH "/Users/michaelrabinovich/ifopt")
#include(ifopt-config)

# Dependencies are linked as INTERFACE targets unless libigl is compiled as a static library
if(LIBIGL_USE_STATIC_LIBRARY)
  set(IGL_SCOPE PUBLIC)
else()
  set(IGL_SCOPE INTERFACE)
endif()

set(LIBIGL_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../libigl/")
set(LIBIGL_SOURCE_DIR "${LIBIGL_ROOT}/include")
set(LIBIGL_EXTERNAL "${LIBIGL_ROOT}/external")


#add_subdirectory(${LIBIGL_EXTERNAL}/libigl-imgui imgui)
#target_link_libraries(igl_opengl_glfw_imgui ${IGL_SCOPE} igl_opengl_glfw imgui)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -ffast-math -O2")

# Add your project files
file(GLOB SRCFILES *.cpp)

file(GLOB SRC_OPTIMIZATION 
      Optimization/*.cpp
      Optimization/Solvers/*.cpp
      Optimization/Solvers/lbfgs/*.cpp)

file(GLOB SRC_QUAD QuadMesh/*.cpp)
file(GLOB SRC_DOG 
    Dog/*.cpp
    Dog/Objectives/*.cpp
    Dog/Solvers/*.cpp)
file(GLOB SRC_CREASE_PATTERNS CreasePatterns/*.cpp)
file(GLOB SRC_GUI Gui/*.cpp)
file(GLOB SRC_GEOMETRY_PRIMITIVES GeometryPrimitives/*.cpp)
file(GLOB SRC_FOLDING Folding/*.cpp)
#file(GLOB IFOPT_SRC ifopt_core/src/*.cc)
#file(GLOB IFOPT_IPOPT_SRC ifopt_ipopt/src/*.cc)

#message("${IFOPT_SRC}")
#find_package(OpenMP)
#find_package(OpenMP3)
if(OpenMP_CXX_FOUND)
    target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
endif()

add_library(libpardiso600-GNU720-X86-64 SHARED IMPORTED)
set_target_properties(libpardiso600-GNU720-X86-64 PROPERTIES IMPORTED_LOCATION "./libpardiso600-GNU720-X86-64.so")

#add_library(libifopt_core SHARED IMPORTED)
#set_target_properties(libifopt_core PROPERTIES IMPORTED_LOCATION "libifopt_core.dylib")

#add_library(libifopt_ipopt SHARED IMPORTED)
#set_target_properties(libifopt_ipopt PROPERTIES IMPORTED_LOCATION "libifopt_ipopt.dylib")

#target_include_directories(ifopt "/Users/michaelrabinovich/ifopt")
#add_definitions("-fopenmp")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")



#add_definitions(-DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
#add_definitions(-DOpenMP_CXX_LIB_NAMES="omp")
#add_definitions(-DOpenMP_omp_LIBRARY=/usr/local/opt/libomp/lib/libomp.dylib)
#add_definitions()
#SET(CMAKE_EXE_LINKER_FLAGS  "-L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include")

#add_library(libomp SHARED IMPORTED)
#set_target_properties(libomp PROPERTIES IMPORTED_LOCATION "libomp.dylib")

add_executable(${PROJECT_NAME}_bin ${SRCFILES} ${SRC_OPTIMIZATION} ${SRC_QUAD} ${SRC_DOG} ${SRC_CREASE_PATTERNS} ${SRC_GUI} ${SRC_GEOMETRY_PRIMITIVES} ${SRC_FOLDING})
target_link_libraries(${PROJECT_NAME}_bin PUBLIC igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::cgal igl::triangle libpardiso600-GNU720-X86-64 gomp gfortran pthread openblas)

#-lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_sequential -lmkl_core -lm
#-lgfortran -lgfortran -lm
