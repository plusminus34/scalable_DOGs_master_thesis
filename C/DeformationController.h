#include "DogEditor.h"


// The DOG structure should be constant (so should have a global DOG)
// Then there should be some pointer to an "active DOG"
// If the active dog is the global one, we just use the same one (everything should be updated on spot)
//	If the active dog is a submesh, we operate on a copy, and allow the user to "merge" these in the GUI
// Every change to this one should create a new DogEditor on that mesh. 
class DeformationController {
public:
	DeformationController(igl::opengl::glfw::Viewer& viewer) : viewer(viewer) {}
	void single_optimization();
	void init_from_new_dog(igl::opengl::glfw::Viewer& viewer, Dog& dog, const QuadTopology& quadTop);

	int edited_mesh = 0; // 0 means the entire mesh, i means the i-1 connected component submesh	

private:
	void reset_editor();
};