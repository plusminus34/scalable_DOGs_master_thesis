#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include "../QuadMesh/Quad.h"
#include "../Dog/Dog.h"

void get_wireframe_edges(const Eigen::MatrixXd& V, const QuadTopology& quadTop, Eigen::MatrixXd& E1, Eigen::MatrixXd& E2, bool display_border);
void render_wireframe(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const QuadTopology& quadTop, bool display_border = true);
void render_dog_stitching_curves(igl::opengl::glfw::Viewer& viewer, const Dog& dog,
					Eigen::RowVector3d color = Eigen::RowVector3d(0, 0, 0));
void render_dog_stitching_constraints(igl::opengl::glfw::Viewer& viewer, const Dog& dog, 
					Eigen::RowVector3d color = Eigen::RowVector3d(0, 0, 0));
void render_wireframe_boundary(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const QuadTopology& quadTop,
		Eigen::RowVector3d color);