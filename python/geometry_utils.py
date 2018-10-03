from shapely.geometry import *
from shapely.geometry.polygon import *
from shapely.prepared import prep
from shapely.ops import cascaded_union,linemerge,split,snap
import numpy as np

def split_line_to_geometry(polygon, line):
	#union = cascaded_union(split(line, polygon))
	split_res = split(line, polygon)
	# If there was a split, unite the lines
	if len(list(split_res)) > 1:
		#print 'splitting!'
		return linemerge(cascaded_union(split_res))
	else: # otherwise return the original line
		return line

def extend_polyline_by_epsilon_forwards(polyline, eps):
	f_point, ff_point = np.array(polyline.coords[0]), np.array(polyline.coords[1])
	diff_vec = f_point-ff_point
	direction_vec = diff_vec/np.linalg.norm(diff_vec)
	new_first_coord = ff_point + diff_vec + eps*direction_vec # extend it by an epsilon
	new_coords = [new_first_coord] + polyline.coords[1:]
	return LineString(new_coords)

def extend_polyline_by_epsilon_backwards(polyline, eps):
	b_point, bb_point = np.array(polyline.coords[-1]), np.array(polyline.coords[-2])
	diff_vec = b_point-bb_point
	direction_vec = diff_vec/np.linalg.norm(diff_vec)
	new_last_coord = bb_point + diff_vec + eps*direction_vec # extend it by an epsilon

	new_coords = polyline.coords[0:-1]
	new_coords.append((new_last_coord[0],new_last_coord[1]))
	#print 'new coords = ', new_coords
	return LineString(new_coords)

# extend a polygon forward and backward by an epsilon so that a 'split' operation will not numerically fail
def extend_polyline_by_epsilon(polyline, eps):
	forward_eps = extend_polyline_by_epsilon_forwards(polyline,eps)
	return extend_polyline_by_epsilon_backwards(forward_eps,eps)
