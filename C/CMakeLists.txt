cmake_minimum_required(VERSION 3.1)
project(dog_editor)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_CGAL            "Use CGAL"  ON)

set(CMAKE_BUILD_TYPE "Release")

find_package(LIBIGL REQUIRED QUIET)

# Add your project files
file(GLOB SRCFILES *.cpp)
file(GLOB SRC_CREASE_PATTERNS CreasePatterns/*.cpp)

add_executable(${PROJECT_NAME}_bin ${SRCFILES} ${SRC_CREASE_PATTERNS})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::cgal)