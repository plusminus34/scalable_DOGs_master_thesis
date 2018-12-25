#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include "Lasso.h"

class Editor {
public:
	Editor(igl::opengl::glfw::Viewer& viewer,
		  const Eigen::MatrixXd &V,
         const Eigen::MatrixXi &F_tri);

	enum MouseMode { SELECT, TRANSLATE, NONE};//, ROTATE, CUT, GLUE1, GLUE2, NONE };
	enum SelectMode {VertexPicker, PathPicker, CurvePicker};

	MouseMode mouse_mode = NONE;
	SelectMode select_mode = VertexPicker;
	bool callback_mouse_down(int button, int modifier);
private:
	void applySelection();
	void onNewHandleID();
	void compute_handle_centroids();
	void get_new_handle_locations();
	Eigen::Vector3f computeTranslation (int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
	Eigen::Vector4f computeRotation(int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);

	igl::opengl::glfw::Viewer& viewer;
	const Eigen::MatrixXd &V;
    const Eigen::MatrixXi &F;
	Lasso lasso;

	Eigen::Vector3f translation;
	Eigen::VectorXi b; Eigen::VectorXd bc;

	//list of currently selected vertices
	Eigen::VectorXi selected_v;//(0,1);
	//updated positions of handle vertices, #HV x3
	Eigen::MatrixXd handle_vertex_positions;//(0,3);
	//list of all vertices belonging to handles, #HV x1
	Eigen::VectorXi handle_vertices;//(0,1);
	//for saving constrained vertices
	//vertex-to-handle index, #V x1 (-1 if vertex is free)
	Eigen::VectorXi handle_id;//(0,1);
	int current_handle = -1;

	//index of handle being moved
	int moving_handle = -1;

	//centroids of handle regions, #H x1
	Eigen::MatrixXd handle_centroids;//(0,3);
	int down_mouse_x = -1, down_mouse_y = -1;

	bool deforming = false;
};
