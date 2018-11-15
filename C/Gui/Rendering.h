#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include "../QuadMesh/Quad.h"
#include "../Dog/Dog.h"

void get_wireframe_edges(const Eigen::MatrixXd& V, const QuadTopology& quadTop, Eigen::MatrixXd& E1, Eigen::MatrixXd& E2);
void render_wireframe(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const QuadTopology& quadTop);
void render_dog_stitching_curves(igl::opengl::glfw::Viewer& viewer, const Dog& dog);