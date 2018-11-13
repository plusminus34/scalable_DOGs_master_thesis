#pragma once

#include "Dog/Dog.h"

#include "igl/serialize.h"

struct ModelState : public igl::Serializable {
	Dog dog;
	QuadTopology quadTop;

	void InitSerialization() {
      Add(dog,std::string("_dog"));
      Add(quadTop,std::string("_quadTop"));
    }
};