#pragma once

#include "Dog/Dog.h"

#include "igl/serialize.h"

struct ModelState : public igl::Serializable {
	Dog dog;
	QuadTopology quadTop;

	void init_from_mesh(const std::string& mesh_path);
	void load_from_workspace(const std::string& workspace_path) {igl::deserialize(*this,"State",workspace_path);}
	void save_to_workspace(const std::string& workspace_path) {igl::serialize(*this,"State",workspace_path);}

	void InitSerialization() {
      Add(dog,std::string("_dog"));
      Add(quadTop,std::string("_quadTop"));
    }
};