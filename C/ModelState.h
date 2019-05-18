#pragma once

#include "Dog/Dog.h"
#include "DeformationController.h"

#include "igl/serialize.h"

struct CreasesVisualization : public igl::Serializable {
	Eigen::MatrixXd V_arr;
	Eigen::MatrixXi F_arr;
	Eigen::MatrixXd faceColors;
	Eigen::MatrixXd edge_pts1,edge_pts2;
	Eigen::MatrixXd meshE1,meshE2;

	void InitSerialization() {
      Add(V_arr,std::string("_V_arr"));
      Add(F_arr,std::string("_F_arr"));
      Add(faceColors,std::string("_faceColors"));
      Add(edge_pts1,std::string("_edge_pts1"));
      Add(edge_pts2,std::string("_edge_pts2"));
      Add(meshE1,std::string("_meshE1"));
      Add(meshE2,std::string("_meshE2"));
    }
};

struct ModelState : public igl::Serializable {
	Dog dog;
	CreasesVisualization creasesVisualization;
	DeformationController DC;

	void init_from_mesh(const std::string& mesh_path);
	void init_from_planar(int square_h, int square_w);
	void init_from_svg(const std::string& svg_path, int x_res, int y_res);
	void load_from_workspace(const std::string& workspace_path) {igl::deserialize(*this,"State",workspace_path);}
	void save_to_workspace(const std::string& workspace_path) {igl::serialize(*this,"State",workspace_path);}

	void InitSerialization() {
      Add(dog,std::string("_dog"));
      Add(creasesVisualization,std::string("_creasesVisualization"));
      Add(DC,std::string("_DC"));
    }
private:
	void setup_non_creased_dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
};

