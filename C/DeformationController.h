#include "Gui/DogEditor.h"


// The DOG structure should be constant (so should have a global DOG)
// Then there should be some pointer to an "active DOG"
// If the active dog is the global one, we just use the same one (everything should be updated on spot)
//	If the active dog is a submesh, we operate on a copy, and allow the user to "merge" these in the GUI
// Every change to this one should create a new DogEditor on that mesh. 
class DeformationController {
public:
	DeformationController() {globalDog = NULL; editedSubmesh = NULL;}
	void init_from_new_dog(Dog& dog);
	void init_viewer(igl::opengl::glfw::Viewer& viewer_i) {viewer = &viewer_i; dogEditor.init_viewer(viewer_i);}

	// pass it on to the editor
	void single_optimization() {return dogEditor.single_optimization();}
	// if needed: change activeDog and update the editor accordingly
	void update_edited_mesh(int newEditedSumeshI);

private:
	igl::opengl::glfw::Viewer* viewer;
	DogEditor dogEditor;
	Dog* globalDog;
	Dog* editedSubmesh;
	int editedSubmeshI = 0; // 0 means the entire mesh, i means the i-1 connected component submesh	
};