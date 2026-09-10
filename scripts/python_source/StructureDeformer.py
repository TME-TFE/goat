# Description
#------------
# Simple module to deform a structure and to move points etc

import numpy as np
# Rotation around point
def RotatePoint(px, py, ox, oy, theta):
    # Description
    #------------
    # Rotate the vector OP for theta amount of degrees around point O. 
    # Returns the new coordinates of the rotated point P. Theta should 
    # be in degrees

    # Initialize
    dx = px - ox 
    dy = py - oy 
    theta = theta/180*np.pi # convert to radians 

    # Rotate
    newpx = ox + dx*np.cos(theta) - dy*np.sin(theta)
    newpy = oy + dx*np.sin(theta) + dy*np.cos(theta)

    return newpx, newpy

# Translation along a vector
def TranslatePoint(px, py, dx, dy, length):
    # Description
    #------------
    # Translate a point by a length 'length' along the vector dx, dy

    # Initialize
    dn = np.sqrt(dx**2 + dy**2)
    dx = dx/dn 
    dy = dy/dn 

    # Translate
    newpx = px + dx*length 
    newpy = py + dy*length 

    return newpx, newpy

# Add a new point in an existing polygon
def AddPoint(px, py, npx, npy, addbefore, pointindex):
    # Description
    #------------
    # Add a point with coordinates npx, npy in the existing polygon 
    # defined by the points px, py before (addbefore >= 0) or after 
    # (addbefore < 0) the point defined by pointindex 

    # Checks
    if len(px) != len(py):
        raise ValueError("AddPoint: px and py should have equal sizes")
    if (pointindex < 0) or (pointindex >= len(px)):
        raise ValueError("AddPoint: pointindex out of range")
    
    # Add point 
    if addbefore < 0:
        # Add after pointindex 
        px = np.insert(px, pointindex+1, npx)
        py = np.insert(py, pointindex+1, npy)
    else:
        px = np.insert(px, pointindex, npx)
        py = np.insert(py, pointindex, npy)

    # Return
    return px, py
    
# Reset the coordinates of a point to be on the intersection of two lines
# defined by their tangential vectors (not normalized)
def LineIntersections(x11, y11, x12, y12, x21, y21, x22, y22):
    disttol = 1e-12
    macheps = 1e-12
    # Compute
    dx1 = x12 - x11 
    dy1 = y12 - y11 
    dx2 = x22 - x21 
    dy2 = y22 - y21 

    # Some checks
    if ( (abs(dx1) < disttol) and (abs(dy1) < disttol) ): 
        # Print warning
        print('LineIntersections: ' +
                'first line is a point up to distance precision, ' +
                'results may be inaccurate. Consider rescaling the' +
                ' coordinates')
        
    if ( (abs(dx2) < disttol) and (abs(dy2) < disttol) ): 
        # Print warning
        print('LineIntersections: ' +
                'second line is a point up to distance precision, ' +
                'results may be inaccurate. Consider rescaling the' +
                ' coordinates')
        

    # Check determinant
    det = -dy1*dx2 + dx1*dy2
    if ( (abs(det) < macheps) ): 
        # Parallel or collinear lines - need to check collinearity
        
        # Compute normal
        nx = -(y11 - y12)
        ny = (x11 - x12)
        nn = np.sqrt(nx**2 + ny**2)

        # Compute vector between lines
        vx = (x11 - x21)
        vy = (y11 - y21)

        # Compute the distance
        dist = abs( vx*nx/nn + vy*ny/nn )

        # Check 
        if (dist < disttol): 
            # collinear lines, return inf
            x = np.nan
            y = x
            return x, y  
        else:
            # Parallel lines, return nan
            x = np.nan
            y = x 
            
            return x, y 
        

    # Compute intersection
    #=====================
    # Hedge for small dx when computing slope 
    if (abs(dx1) > abs(dx2)): 
        if (abs(dx1) > disttol): 
            r1 = dy1/dx1 
            if (abs(dx2) > disttol): 
                # Two non-parallel, non-vertical and non-horizontal lines
                r2 = dy2/dx2 
                x = (r1*x11 - r2*x21 -y11 + y21)/(r1 - r2)
                y = r1*(x - x11) + y11
            else:
                # Second line is vertical line, first one is non-vertical
                x = x21
                y = r1*(x - x11) + y11
                

        else:
            # Both lines are parallel - should've been captured before actually
            # Compute normal
            nx = -(y11 - y12)
            ny = (x11 - x12)
            nn = np.sqrt(nx**2 + ny**2)

            # Compute vector between lines
            vx = (x11 - x21)
            vy = (y11 - y21)

            # Compute the distance
            dist = abs( vx*nx/nn + vy*ny/nn )

            # Check 
            if (dist < disttol): 
                # collinear lines, return inf
                x = np.nan
                y = x 
                return x, y 
            else:
                # Parallel lines, return nan
                x = np.nan
                y = x 
                
                return x, y
            
    else:
        if (abs(dx2) > disttol): 
            r2 = dy2/dx2
            if (abs(dx1) > disttol): 
                # Two non-parallel, non-vertical and non-horizontal lines
                r1 = dy1/dx1 
                x = (r1*x11 - r2*x21 -y11 + y21)/(r1 - r2)
                y = r2*(x - x21) + y21
            else:
                # First line is vertical line, second one is non-vertical
                x = x11
                y = r2*(x - x21) + y21
                

        else:
            # Both lines are parallel - should've been captured before actually
            # Compute normal
            nx = -(y11 - y12)
            ny = (x11 - x12)
            nn = np.sqrt(nx**2 + ny**2)

            # Compute vector between lines
            vx = (x11 - x21)
            vy = (y11 - y21)

            # Compute the distance
            dist = abs( vx*nx/nn + vy*ny/nn )

            # Check 
            if (dist < disttol): 
                # collinear lines, return inf
                x = np.nan
                y = x 
            else:
                # Parallel lines, return nan
                x = np.nan
                y = x 
                
    return x, y
            
        
