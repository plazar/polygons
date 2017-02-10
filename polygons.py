import shapefile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


import pyximport
pyximport.install()
from contains_point import contains_point as contains_point_cython

class Polygon(object):
    def __init__(self, verts, info, shapefn=None, shapenum=None, polynum=None):
        self.verts = np.array(verts)
        assert self.verts.ndim == 2
        self.info = info
        self.shapefn = shapefn
        self.shapenum = shapenum
        self.polynum = polynum
        self.area = self.__polygon_signed_area()
        self.centroid = self.__polygon_centroid()
        
    def __polygon_signed_area(self):
        # The signed area of a polygon is positive if
        # vertices are ordered in CCW order, and
        # negative is ordered in CW order
        x, y = np.array(self.verts).T
        x1 = np.roll(x, -1)
        y1 = np.roll(y, -1)
        return 0.5*np.sum(x*y1 - x1*y)

    def __polygon_centroid(self):
        x, y = np.array(self.verts).T
        x1 = np.roll(x, -1)
        y1 = np.roll(y, -1)
        A = self.__polygon_signed_area()
        centx = np.sum((x + x1)*(x*y1 - x1*y))/(6.0*A)
        centy = np.sum((y + y1)*(x*y1 - x1*y))/(6.0*A)
        return centx, centy

    def contains_points(self, xys, *args, **kwargs):
        """Check if a points are within the polygon.

            Inputs:
                xys: A 2D array of (x,y) coordinates for each points
                ** Additional positional and keyword arguments are
                    passed on to "self.contains_point(...)"
            Output:
                is_contained: A numpy array of boolean values -
                    True if the point is contained, False otherwise.
        """
        return np.apply_along_axis(self.contains_point, axis=1, arr=xys, 
                                   *args, **kwargs)

    def contains_point(self, *args, **kwargs):
        return contains_point_cython(self.verts, *args, **kwargs)

    def contains_point_old(self, xy, vertex_is_in=True, edge_is_in=True):
        """Check if a point is within the polygon.
            This code is adapted from: 
                http://geospatialpython.com/2011/08/point-in-polygon-2-on-line.html

            Inputs:
                xy: A tuple of (x,y) coordinates of the point
                vertex_is_in: If True, point on a vertex is considered to be contained.
                    (Default: True)
                edge_is_in: If True, point on an edge is considered to be contained.
                    (Default: True)

            Output:
                is_contained: True if the point is contained, False otherwise.
        """
        x, y = xy
        # check if point is a vertex
        if np.any(np.bitwise_and((x == self.verts[:,0]), (y == self.verts[:,1]))): 
            return vertex_is_in
 
        nverts = len(self.verts)
        # check if point is on a boundary
        for ii in xrange(nverts):
            p1 = None
            p2 = None
            if ii == 0:
                p1 = self.verts[0]
                p2 = self.verts[1]
            else:
                p1 = self.verts[ii-1]
                p2 = self.verts[ii]
            if (p1[1] == p2[1]) and (p1[1] == y) and \
                    (x > min(p1[0], p2[0])) and (x < max(p1[0], p2[0])):
                return edge_is_in
           
        inside = False 
        p1x, p1y = self.verts[0]
        for ii in xrange(nverts+1):
            p2x, p2y = self.verts[ii % nverts]
            if min(p1y, p2y) < y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
            p1x, p1y = p2x, p2y
        return inside


def plot_polygons(polygons, constraint=None, label_index=None, 
                  label_size='small', label_colour='r', 
                  maxpoly=None, indices=None, **kwargs):
    """Plot a shapefile as polygons to the current matplotlib axes.
    
        Inputs:
            polygons: A list of Polygon objects.
            constraint: a function that when provided with the
                        Polygon object will return a boolean value 
                        specifying whether or not the shape should 
                        be plotted (Default: show all)
            label_index: the index of the Polygon's info list to use 
                         as a label. If None, don't show labels. 
                         If 'polynum', show polynomial number.
                         If 'shapenum', show shape number.
                         (Default: None - i.e. no label)
            label_size: Size of labels (Default: 'small')
            maxpoly: The maximum number of polynomials to show.
            indices: The indices of the shapes to show.
            
            ** additional keyword arguments are passed directly to the 
               matplotlib Polygon patches **
    """
    ax = plt.gca() # Get current axes or create one
    polystyle = dict(lw='2', fc='none', ec='k')
    polystyle.update(kwargs)
 
    if constraint is None:
        constraint = lambda poly: True
    for poly in polygons:
        if (maxpoly is not None) and (npoly > maxpoly):
            break
        if not constraint(poly):
            continue
        polypatch = plt.Polygon(poly.verts, **polystyle)
        ax.add_patch(polypatch)
        centx, centy = poly.centroid
        if label_index is not None:
            if label_index == 'polynum':
                lbl = str(poly.polynum)
            elif label_index == 'shapenum':
                lbl = str(poly.shapenum)
            else:
                lbl = poly.info[label_index]
            plt.text(centx, centy, lbl, fontsize=label_size, clip_on=True,
                     ha='center', va='center', fontdict=dict(color=label_colour))
        if (maxpoly is not None) and (npoly > maxpoly):
            break


def polygons_from_shapefile(shapefn, constraint=None, maxpoly=None):
    """Get polygons from the shapefile provided.
        
        Inputs:
            shapefn: The name of the shapefile.
            constraint: a function that when provided with the
                        Polygon object will return a boolean value 
                        specifying whether or not the shape should 
                        be kept. (Default: keep all)
            maxpoly: The maximum number of polynomials to keep.

        Outputs:
            polys: A list of Polygon objects
    """

    sf = shapefile.Reader(shapefn)
    polygons = []
    
    if constraint is None:
        constraint = lambda poly: True

    npoly = 0
    nshape = 0

    indices = xrange(sf.numRecords)
    for ii in indices:
        shp = sf.shape(ii)
        rec = sf.record(ii)
        if not hasattr(shp, 'parts'):
            continue
        inds = list(shp.parts)+[-1]
        start = inds[0]
        for stop in inds[1:]:
            pts = shp.points[start:stop]
            poly = Polygon(pts, rec, shapefn=shapefn, 
                           shapenum=nshape, polynum=npoly)
            # Check if polygon should be kept
            if (maxpoly is not None) and (npoly > maxpoly):
                break
            if not constraint(poly):
                continue
            polygons.append(poly)
            start = stop
            npoly += 1
        nshape += 1
    return polygons
