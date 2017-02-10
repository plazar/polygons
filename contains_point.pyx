def contains_point(verts, xy, vertex_is_in=True, edge_is_in=True):
    """Check if a point is within the polygon defined by verts.
        This code is adapted from: 
            http://geospatialpython.com/2011/08/point-in-polygon-2-on-line.html

        Inputs:
            verts: a list of vertices
            xy: A tuple of (x,y) coordinates of the point
            vertex_is_in: If True, point on a vertex is considered to be contained.
                (Default: True)
            edge_is_in: If True, point on an edge is considered to be contained.
                (Default: True)

        Output:
            is_contained: True if the point is contained, False otherwise.
    """
    x, y = xy
    
    v1x = verts[-1][0]
    v1y = verts[-1][1]
    inside = False
    for vert in verts:
        if (x == vert[0]) and (y == vert[1]):
            inside = vertex_is_in
            break
        v2x = vert[0]
        v2y = vert[1]
        
        if v1y > v2y:
            ymax = v1y
            ymin = v2y
        else:
            ymax = v2y
            ymin = v1y

        if v1x > v2x:
            xmax = v1x
            xmin = v2x
        else:
            xmax = v2x
            xmin = v1x
        # Check if point is on the boundary
        if (v2x == v1x):
            if (y > ymin) and (y < ymax):
                inside = edge_is_in
                break
            else:
                continue
        m = (v2y-v1y)/(v2x-v1x)
        b = v2y-m*v2x
        yedge_at_x = (m*x+b)
        if (x > xmin) and (x < xmax):
            if (y == yedge_at_x):
                inside = edge_is_in
                break
            elif (y < yedge_at_x):
                inside = not inside
        # Update point
        v1x = v2x
        v1y = v2y
    return inside

