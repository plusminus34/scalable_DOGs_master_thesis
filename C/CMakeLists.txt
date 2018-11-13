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

### Compile the ImGui part ###
#if(LIBIGL_WITH_OPENGL_GLFW_IMGUI)
#  if(TARGET igl::opengl_glfw)
    # ImGui module
    compile_igl_module("opengl/glfw/imgui")
    if(NOT TARGET imgui)
      igl_download_imgui()
      add_subdirectory(${LIBIGL_EXTERNAL}/libigl-imgui imgui)
    endif()
    target_link_libraries(igl_opengl_glfw_imgui ${IGL_SCOPE} igl_opengl_glfw imgui)
 # endif()
#endif()

#add_subdirectory(${LIBIGL_EXTERNAL}/libigl-imgui imgui)
#target_link_libraries(igl_opengl_glfw_imgui ${IGL_SCOPE} igl_opengl_glfw imgui)

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

add_library(libpardiso500-MACOS-X86-64 SHARED IMPORTED)
set_target_properties(libpardiso500-MACOS-X86-64 PROPERTIES IMPORTED_LOCATION "libpardiso500-MACOS-X86-64.dylib")
add_executable(${PROJECT_NAME}_bin ${SRCFILES} ${SRC_OPTIMIZATION} ${SRC_QUAD} ${SRC_DOG} ${SRC_CREASE_PATTERNS} ${SRC_GUI})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::cgal igl::triangle libpardiso500-MACOS-X86-64)