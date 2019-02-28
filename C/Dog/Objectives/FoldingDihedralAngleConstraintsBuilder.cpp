#include "FoldingDihedralAngleConstraintsBuilder.h"

FoldingDihedralAngleConstraintsBuilder::FoldingDihedralAngleConstraintsBuilder(const Dog& dog, 
			std::vector<double>& destination_dihedral_angles, const double& timestep) : 
				dog(dog), destination_dihedral_angles(destination_dihedral_angles), timestep(timestep) {
	// Todo: find the initial tangent angles so that we could convert angles to angles
}