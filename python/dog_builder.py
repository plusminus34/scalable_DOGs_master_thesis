import numpy as np
import os, sys
sys.path.insert(0, os.getcwd() + "/../../libigl/python")
sys.path.insert(0, os.getcwd() + "/../../libigl/external/nanogui/build/python")
import math
import pyigl as igl

from polygons_to_orthogonal_grids import *

def dog_builder(V_list,F_list, polylines):
	# return global V,F, then V and F for rendering, and constraints 

	# The global F just concatenates the vertices and faces list (after adding indices to faces)
	#
	# The constraints data is defined per vertex on the polylines.
	# A polyline vertex is defined uniquely by its 2D coordinates (so there won't be doubles). 
	# No need to save it but this is a good way to filter multiples of those.
	# Each one of those is either exactly a vertex, or a point on an edge.
	# In any case it can be represented as a point of an edge. So for every polygon that does contain this point:
	#	We save two indices of vertices and an appropriate weight
	#   If we have 2 pairs of this we need one constraint, if we have 3 pairs 2, and in general for n we need n-1 (1 to 2, 2 to 3, etc)
	#   So at the end we just need to save a list of constraints of the form:
	#				2 points on edges being equal: (v_i1,v_i2) 0<=w_1<=1, (v_j1,v_j2) 0<=w_2<=1

	# So the function should find all unique points lying on polylines, then find all polygons that touch those,
	# 	then find these as points on (some) edges for each one (n such edges), and save those
	#										 (the optimization code could use them as n-1 constraints)
	#
	#
	# Rendering requires culling the faces. 
	# The triangulation can be an arbitrary one (triangulate) if you take the unique vertices locations together with the polylines.
	# If we have a unique order for the polylines vertices (given by another list of edges such that no vertex is there twice)
	#   then we just need the list of edges.
	# So basically: Save a list of vertices indices in the global mesh (for the inner vertices positions)
	# 				A list of edges for polylines (for the polylines)
	#				And an F for triangulation that assumes this order of vertices (so for initialization build it and triangulate)
	#
	# In case the triangulation looks bad, do a manual one
	V,F = get_global_VF(V_list, F_list)
	pass

def get_global_VF(V_list, F_list):
	V,F = V_list[0], F_list[0]
	max_f = F.max()
	for mesh_i in range(1, len(V_list)):
		V = np.concatenate((V, V_list[mesh_i]))
		F = np.concatenate((F, F_list[mesh_i] + F.max() + 1))
	return V,F

def test_dog_builder(svg_file):
	border_poly,polylines = svg_creases_to_polygonal_data(svg_file)
	face_polygons, polylines = crease_pattern(border_poly, polylines)

	res_x, res_y = 10,10
	V_list, F_list, polylines = polygons_to_orthogonal_grids(face_polygons, border_poly, polylines, res_x, res_y)
	dog_builder(V_list, F_list, polylines)

if __name__ == "__main__":
	if len(sys.argv) > 1:
		test_dog_builder(sys.argv[1])
	else:
		test_dog_builder("../crease_patterns/1_curve.svg")
