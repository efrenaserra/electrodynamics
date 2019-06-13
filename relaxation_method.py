# -*- coding: utf-8 -*-
"""
Solution to Problem 3.30 of Chapter 3, Electric Fields Around Conductors
in Electricity and Magnetism by Edward M. Purcell, 2nd Edition.

@author: Efren A. Serra
"""
import numpy as np
from enum import Enum
import pylab

""" Boundary or interior point """
class point_t(Enum):
    boundary = 0
    interior = 1

""" Grid point """
class grid_point(object):
    # Constructor
    def __init__(self, point_t=point_t.interior):
        self.point_t = point_t

    def __repr__(self):
        return "<grid_point: point_t:%s>" % self.point_t


""" Rectangular grid point """
class rect_grid_point(grid_point):
    # Constructor
    def __init__(self, i, j, data):
        super().__init__()
        self.i = i
        self.j = j
        self.data = data

    def __repr__(self):
        return "<rect_grid_point: i:%d j:%d data:%f point_t:%s>" % (self.i, self.j, self.data, self.point_t)

    def set_point_t(self, point_t):
        self.point_t = point_t

    # Relax
    def relax(self, nearest_n):
        relaxed_val = self.data
        for n in nearest_n:
            relaxed_val += n.data()

        relaxed_val /= len(nearest_n)
        return relaxed_val

""" Rectangular grid """
class rect_grid(object):
    # Constructor
    def __init__(self, nx, ny, inner_x, inner_nx, inner_y, inner_ny):
        self.nx       = nx
        self.ny       = ny
        self.inner_x  = inner_x
        self.inner_nx = inner_nx
        self.inner_y  = inner_y
        self.inner_ny = inner_ny
        self.xy_point = np.array([[rect_grid_point(i, j, 0.) for j in range(ny)] for i in range(nx)])

    def __repr__(self):
        data = np.array([[self.xy_point[i,j].data for j in range(self.ny)] for i in range(self.nx)], dtype=np.float32)
        return "<rect_grid: nx: %d ny: %d data:\n%s>" % (self.nx, self.ny, data)

    def data(self):
        data = np.array([[self.xy_point[i,j].data for j in range(self.ny)] for i in range(self.nx)], dtype=np.float32)
        return data

    # Sets outer boundary grid values
    def set_outer(self, val):
        for bp in self.xy_point[0,:]:
            bp.data = val
            bp.set_point_t(point_t.boundary)
        # Top row
        for bp in self.xy_point[self.nx-1,:]:
            bp.data = val
            bp.set_point_t(point_t.boundary)
        # Bottom row
        for bp in self.xy_point[:,0]:
            bp.data = val
            bp.set_point_t(point_t.boundary)
        # Left column
        for bp in self.xy_point[:,self.ny-1]:
            bp.data = val
            bp.set_point_t(point_t.boundary)
        # Right column

    # Sets inner boundary grid values
    def set_inner(self, val):
        x = self.inner_x
        y = self.inner_y
        inner_nx = self.inner_nx
        inner_ny = self.inner_ny

        for bp in self.xy_point[x,y:y+inner_ny-1]:
            bp.data = val
            bp.set_point_t(point_t.boundary)
        # Top

        for bp in self.xy_point[x+inner_nx-1,y:y+inner_ny-1]:
            bp.data = val
            bp.set_point_t(point_t.boundary)
        # Bottom

        for bp in self.xy_point[x:x+inner_nx,y]:
            bp.data = val
            bp.set_point_t(point_t.boundary)
        # Left

        for bp in self.xy_point[x:x+inner_nx,y+inner_ny-1]:
            bp.data = val
            bp.set_point_t(point_t.boundary)
        # Right

    def relax(self):
        for j in range(self.ny):
            for i in range(self.nx):
                if self.xy_point[i,j].point_t == point_t.boundary:
                    pass
                elif (self.inner_x <= i and i <= (self.inner_x + self.inner_nx - 1)) and \
                (self.inner_y <= j and j <= (self.inner_y + self.inner_nx - 1)):
                    pass
                elif self.xy_point[i,j].point_t == point_t.interior:
                    self.xy_point[i,j].data =     \
                    (self.xy_point[i-1, j].data + \
                     self.xy_point[i, j+1].data + \
                     self.xy_point[i, j-1].data + \
                     self.xy_point[i+1, j].data)/4.

        data = np.array([[self.xy_point[i,j].data for j in range(self.ny)] for i in range(self.nx)], dtype=np.float32)
        return sum(data)

            
def relaxation_method():
    rg = rect_grid(10, 10, 3, 4, 3, 4)

    rg.set_outer(0.0)
    rg.set_inner(100.0)

    n = 0
    n_runs = 10**3
    while (n <= n_runs):
        print (rg)
        print (rg.relax())
        n += 1

    pylab.imshow(rg.data(), origin='lower', cmap='RdGy')
    pylab.colorbar()

def main():
    np.set_printoptions(precision=1)
    relaxation_method()

if __name__ ==  "__main__":
    main()