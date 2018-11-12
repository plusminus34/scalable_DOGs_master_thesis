cmake_minimum_required(VERSION 3.1)
project(dog_editor)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_CGAL            "Use CGAL"  ON)
option(LIBIGL_WITH_TRIANGLE            "Use Triangle"  ON)

set(CMAKE_BUILD_TYPE "Release")

find_package(LIBIGL REQUIRED QUIET)

# Dependencies are linked as INTERFACE targets unless libigl is compiled as a static library
if(LIBIGL_USE_STATIC_LIBRARY)
  set(IGL_SCOPE PUBLIC)
else()
  set(IGL_SCOPE INTERFACE)
endif()

set(LIBIGL_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../libigl/")
set(LIBIGL_SOURCE_DIR "${LIBIGL_ROOT}/include")
set(LIBIGL_EXTERNAL "${LIBIGL_ROOT}/external")

#if(LIBIGL_WITH_TRIANGLE)
  set(TRIANGLE_DIR "${LIBIGL_EXTERNAL}/triangle")
  if(NOT TARGET triangle)
    igl_download_triangle()
    add_subdirectory("${TRIANGLE_DIR}" "triangle")
  endif()
  compile_igl_module("triangle")
  target_link_libraries(igl_triangle ${IGL_SCOPE} triangle)
  target_include_directories(igl_triangle ${IGL_SCOPE} ${TRIANGLE_DIR})
#endif()

# Add your project files
file(GLOB SRCFILES *.cpp)
file(GLOB SRC_OPTIMIZATION Optimization/*.cpp)
file(GLOB SRC_QUAD QuadMesh/*.cpp)
file(GLOB SRC_DOG Dog/*.cpp)
file(GLOB SRC_DOG_OBJECTIVES Dog/Objectives/*.cpp)
file(GLOB SRC_CREASE_PATTERNS CreasePatterns/*.cpp)

add_executable(${PROJECT_NAME}_bin ${SRCFILES} ${SRC_OPTIMIZATION} ${SRC_QUAD} ${SRC_DOG} ${SRC_DOG_OBJECTIVES} ${SRC_CREASE_PATTERNS})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::cgal igl::triangle)