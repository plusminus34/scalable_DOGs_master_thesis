from matplotlib import pyplot
from shapely.geometry import LineString
from descartes import PolygonPatch
import numpy as np

COLOR = {
    True:  '#6699cc',
    False: '#ffcc33'
    }

def v_color(ob):
    return COLOR[ob.is_simple]

def plot_coords(ax, ob, col = '#999999'):
    x, y = ob.xy
    ax.plot(x, y, 'o', color=col, zorder=1)

def plot_bounds(ax, ob):
    x, y = zip(*list((p.x, p.y) for p in ob.boundary))
    ax.plot(x, y, 'o', color='#000000', zorder=1)

def plot_line(ax, ob, l_width = 3, col=''):
    x, y = ob.xy
    if col == '':
        col = v_color(ob)
    ax.plot(x, y, color=col, alpha=0.7, linewidth=l_width, solid_capstyle='round', zorder=2)

def plot_lines(ax, lines):
    for l in lines:
        plot_line(ax,l)

def plot_border_polygon_and_lines(ax, border_polygon, polylines, title = ''):
    # plot border polygon
    pol_patch = PolygonPatch(border_polygon)
    ax.add_patch(pol_patch)
    ax.set_title(title)
    # plot lines
    for line in polylines:
        plot_line(ax, line)
        plot_coords(ax, line)

def plot_face_polygons(face_polygons, polylines, ax, title = ''):
    pol_colors = get_spaced_colors(len(face_polygons))
    i = 0
    for pol in face_polygons:
        color = np.array(pol_colors[i])/255.
        #print 'color = ', color
        #print 'pol.area = ', pol.area
        pol_patch = PolygonPatch(pol, facecolor=color)
        #print 'pol_patch = ', pol_patch
        ax.add_patch(pol_patch)
        i += 1
    
    """
    # plot lines
    for line in polylines:
        plot_line(ax, line)
        #plot_coords(ax, line)
        pass
    """
    ax.set_title(title)

def plot_grid(grid, ax, l_width = 3, color = ''):
    for l in grid:
        plot_line(ax, l, l_width, color)
        #plot_coords(ax, l)

def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]
    
    return [(int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]