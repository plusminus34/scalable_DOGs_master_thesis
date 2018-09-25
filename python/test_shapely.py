import numpy as np
from shapely.geometry import *
from shapely.geometry.polygon import *
from pprint import pprint
from shapely.ops import polygonize,polygonize_full,linemerge

from matplotlib import pyplot as plt
import matplotlib
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch

def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
    
    return [(int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]

if __name__ == "__main__":
	lines = [
		((0, 0), (0, 1)),
		((0, 0), (1, 0)),
		#((0, 0), (1, 1)),
		((1, 0), (0, 1)),
		((1, 0), (1, 1)),
		((0, 1), (1, 1)),

		]

	new_lines = linemerge(lines)
	for line in new_lines:
		print 'line with ', line
	pprint(list(linemerge(lines)))
	pol_list = list(polygonize(lines))
	result, dangles, cuts, invalids = polygonize_full(lines)
	print 'cuts = ', cuts
	print 'dangles = ', dangles
	print 'invalids = ', invalids
	pprint(pol_list)
	fig = plt.figure(1, figsize=(5,5), dpi=90)
	ax = fig.add_subplot(111)
	pol_colors = get_spaced_colors(len(pol_list))
	i = 0
	for pol in pol_list:
		color = np.array(pol_colors[i])/255.
		print 'color = ', color
		pol_patch = PolygonPatch(pol, facecolor=color)
		ax.add_patch(pol_patch)
		i += 1
	plt.show()
