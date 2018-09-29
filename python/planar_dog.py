from shapely.geometry import *
from shapely.geometry.polygon import *

def grid_from_boundary(bounds, res_x = 20, res_y = 20):
	minx, miny, maxx, maxy = bounds
	step_x = (maxx-minx)/res_x
	step_y = (maxy-miny)/res_y
	grid_lines = []
	for y_i in range(res_y+1):
		points = [(minx + step_x*x_i,miny+y_i*step_y) for x_i in range(res_x+1)]
		grid_lines.append(LineString(points))

	for x_i in range(res_x+1):
		points = [(minx + step_x*x_i,miny+y_i*step_y) for y_i in range(res_y+1)]
		grid_lines.append(LineString(points))

	return grid_lines
