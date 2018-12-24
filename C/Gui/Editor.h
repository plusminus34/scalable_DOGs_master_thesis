#pragma once

//#include "ModelState.h"
//#include "DeformationController.h"
#include <igl/opengl/glfw/Viewer.h>

class Editor {
public:
	//Editor(const Dog& dog);

	enum MouseMode { SELECT, TRANSLATE, NONE};//, ROTATE, CUT, GLUE1, GLUE2, NONE };
	enum SelectMode {VertexPicker, PathPicker, CurvePicker};

	MouseMode mouse_mode = NONE;
	SelectMode select_mode = VertexPicker;
private:

	Eigen::Vector3f computeTranslation (igl::opengl::glfw::Viewer& viewer,
																		int mouse_x,
																		int from_x,
																		int mouse_y,
																		int from_y,
																		Eigen::RowVector3d pt3D);
	Eigen::Vector4f computeRotation(igl::opengl::glfw::Viewer& viewer,
																int mouse_x,
																int from_x,
																int mouse_y,
																int from_y,
																Eigen::RowVector3d pt3D);

	//list of currently selected vertices
	Eigen::VectorXi selected_v;//(0,1);
	//updated positions of handle vertices, #HV x3
	Eigen::MatrixXd handle_vertex_positions;//(0,3);
	//list of all vertices belonging to handles, #HV x1
	Eigen::VectorXi handle_vertices;//(0,1);
	//for saving constrained vertices
	//vertex-to-handle index, #V x1 (-1 if vertex is free)
	Eigen::VectorXi handle_id;//(0,1);

	//index of handle being moved
	int moving_handle = -1;

	//centroids of handle regions, #H x1
	Eigen::MatrixXd handle_centroids;//(0,3);

	//const Dog dog;
	//const DeformationController& deformationController;
};
