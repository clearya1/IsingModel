#!/usr/bin/python
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pylab as plt
import matplotlib.animation as animation
import numpy as np
#from scipy.constants import k
import math
import random
from mpl_toolkits import mplot3d

class IsingPlotter:
    """A class for plotting the Ising Model"""

    def __init__(self, rows, temp, h, shape = "Square", jtype = 1, bc = True, n = "nn"):
        """Initialises the Ising Model"""
        self.rows = rows
        self.cols = rows
        self.temp = temp
        self.shape = shape
        if self.shape == "Square":
            self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols))
            self.dim = self.rows*self.cols
            self.position = np.array([[a,b] for a in np.linspace(0,self.rows-1,self.rows) for b in np.linspace(0,self.cols-1,self.cols)])
        
        if self.shape == "Tri":
            self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols))
            self.dim = self.rows*self.cols
            self.position = np.array([[a * np.sqrt(3)/2,b + (a%2) * 0.5] for a in np.linspace(0,self.rows-1,self.rows) for b in np.linspace(0,self.cols-1,self.cols)])
        
        if self.shape == "3D":
            self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols, self.rows))
            self.dim = self.rows*self.cols*self.rows
            self.position = np.array([[a,b,c] for a in np.linspace(0,self.rows-1,self.rows) for b in np.linspace(0,self.cols-1,self.cols) for c in np.linspace(0,self.rows-1,self.rows)])
        
        if self.shape == "Random0" or self.shape == "Random":
            if self.shape == "Random0":
                self.grid = np.random.choice((-1, 0, 1), size=(self.rows, self.cols))
                self.dim = np.sum(np.abs(self.grid))
            if self.shape == "Random":
                self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols))
                self.dim = np.sum(np.abs(self.grid))
            self.disordergrid = []
            self.position = np.array([[a,b] for a in np.linspace(0,self.rows-1,self.rows) for b in np.linspace(0,self.cols-1,self.cols)])
            for i in range(self.rows):
                for j in range(self.cols):
                    x = random.random()*random.choice([-1,1])
                    y = random.random()*random.choice([-1,1])
                    self.disordergrid.append([x,y])
                    self.position[i*self.rows+j][0]+=x
                    self.position[i*self.rows+j][1]+=y
        
        if self.shape == "Random3D0" or self.shape == "Random3D":
            if self.shape == "Random3D0":
                self.grid = np.random.choice((-1, 0, 1), size=(self.rows, self.cols, self.rows))
                self.dim = np.sum(np.abs(self.grid))
            if self.shape == "Random3D":
                self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols, self.rows))
                self.dim = np.sum(np.abs(self.grid))
            self.disordergrid = []
            self.position = np.array([[a,b,c] for a in np.linspace(0,self.rows-1,self.rows) for b in np.linspace(0,self.cols-1,self.cols) for c in np.linspace(0,self.rows-1,self.rows)])
            for i in range(self.rows):
                for j in range(self.cols):
                    for k in range(self.rows):
                        x = random.random()*random.choice([-1,1])
                        y = random.random()*random.choice([-1,1])
                        z = random.random()*random.choice([-1,1])
                        self.disordergrid.append([x,y,z])
                        self.position[i*self.rows*self.cols+j*self.rows+k][0]+=x
                        self.position[i*self.rows*self.cols+j*self.rows+k][1]+=y
                        self.position[i*self.rows*self.cols+j*self.rows+k][2]+=z

        if self.shape == "NiO":
            self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols))
            for i in range(self.rows):
                for j in range(self.cols):
                    if (i%2 == 0 and j%2 == 1) or (i%2 == 1 and j%2 == 0):
                        self.grid[i,j] = 0
                    self.dim = self.rows*self.cols
                    self.position = np.array([[a,b] for a in np.linspace(0,self.rows-1,self.rows) for b in np.linspace(0,self.cols-1,self.cols)])

        self.counter = 0
        self.counter_list = []
        self.bc = bc
        self.n = n
        self.h = h
        self.jtype = jtype
        
        if self.jtype == 1: #Default set-up for jgrid
            if self.shape == "3D":
                self.jgrid = np.ones([self.rows, self.cols, self.rows, self.rows, self.cols, self.rows])
            else:
                self.jgrid = np.ones([self.rows, self.cols, self.rows, self.cols])
        
        if self.jtype == -1:
            self.jgrid = -1*np.ones([self.rows, self.cols, self.rows, self.cols])

        if self.jtype == 2: #Disordered Coupling
            self.jgrid = np.ones([self.rows, self.cols, self.rows, self.cols])
            for i in range(self.rows):
                for j in range(self.cols):
                    for x in range(self.rows):
                        for y in range(self.cols):
                            self.jgrid[i, j, x, y] = random.random()
                            self.jgrid[x, y, i, j] = self.jgrid[i, j, x, y]

        if self.jtype == -2: #Negative Disordered Coupling
            self.jgrid = np.ones([self.rows, self.cols, self.rows, self.cols])
            for i in range(self.rows):
                for j in range(self.cols):
                    for x in range(self.rows):
                        for y in range(self.cols):
                            self.jgrid[i, j, x, y] = -1*random.random()
                            self.jgrid[x, y, i, j] = self.jgrid[i, j, x, y]
    
        if self.jtype == 3: #Disordered Coupling
            self.jgrid = np.ones([self.rows, self.cols, self.rows, self.cols])
            for i in range(self.rows):
                for j in range(self.cols):
                    for x in range(self.rows):
                        for y in range(self.cols):
                            self.jgrid[i, j, x, y] = random.random()*random.choice(-1,1)
                            self.jgrid[x, y, i, j] = self.jgrid[i, j, x, y]

        if self.jtype == 4: #J = exp(-r)
            self.jgrid = np.ones([self.rows, self.cols, self.rows, self.cols])
            for i in range(self.rows):
                for j in range(self.cols):
                    for x in range(self.rows):
                        for y in range(self.cols):
                            self.jgrid[i, j, x, y] = np.exp(-np.sqrt((x-i)**2+(y-j)**2))

        if self.shape == "Random" or self.shape == "Random0": #J = exp(-r)
            self.jgrid = np.ones([self.rows, self.cols, self.rows, self.cols])
            for i in range(self.rows):
                for j in range(self.cols):
                    for x in range(self.rows):
                        for y in range(self.cols):
                            xr = self.disordergrid[x*self.rows+y][0]
                            yr = self.disordergrid[x*self.rows+y][1]
                            ir = self.disordergrid[i*self.rows+j][0]
                            jr = self.disordergrid[i*self.rows+j][1]
                            self.jgrid[i, j, x, y] = np.exp(-np.sqrt((x+xr-i-ir)**2+(y+yr-j-jr)**2))

        if self.shape == "Random3D" or self.shape == "Random3D0": #J = exp(-r)
            self.jgrid = np.ones([self.rows, self.cols, self.rows, self.rows, self.cols, self.rows])
            for i in range(self.rows):
                for j in range(self.cols):
                    for k in range(self.rows):
                        for x in range(self.rows):
                            for y in range(self.cols):
                                for z in range(self.rows):
                                    xr = self.disordergrid[x*self.rows*self.cols+y*self.rows+z][0]
                                    yr = self.disordergrid[x*self.rows*self.cols+y*self.rows+z][1]
                                    zr = self.disordergrid[x*self.rows*self.cols+y*self.rows+z][2]
                                    ir = self.disordergrid[i*self.rows*self.cols+j*self.rows+k][0]
                                    jr = self.disordergrid[i*self.rows*self.cols+j*self.rows+k][1]
                                    kr = self.disordergrid[i*self.rows*self.cols+j*self.rows+k][2]
                                    self.jgrid[i, j, k, x, y, z] = np.exp(-np.sqrt((x+xr-i-ir)**2+(y+yr-j-jr)**2+(z+zr-k-kr)**2))
        
        if self.n == "fn" and self.bc == True:
            list = [1,2,3]
            self.neighgrid = []
            for i in range(self.rows):
                for j in range(self.cols):
                    left = (i, j-random.choice(list))
                    right = (i, (j+random.choice(list))%self.cols)
                    up = (i-random.choice(list),j)
                    down = ((i+random.choice(list))%self.rows, j)
                    self.neighgrid.append([left, right, up, down])

        if self.shape == "Random" or self.shape == "Random0":
            self.neighgrid = []
            for i in range(self.rows):
                for j in range(self.cols):
                    neighlist = []
                    radiuslist = []
                    for x in range(self.rows):
                        for y in range(self.cols):
                            if self.grid[x,y] != 0 and (i != x or j != y):
                                xr = self.disordergrid[x*self.rows+y][0]
                                yr = self.disordergrid[x*self.rows+y][1]
                                ir = self.disordergrid[i*self.rows+j][0]
                                jr = self.disordergrid[i*self.rows+j][1]
                                radiuslist.append(np.sqrt((x+xr-i-ir)**2+(y+yr-j-jr)**2))
                    radiuslist.sort()
                    for x in range(self.rows):
                        for y in range(self.cols):
                            xr = self.disordergrid[x*self.rows+y][0]
                            yr = self.disordergrid[x*self.rows+y][1]
                            ir = self.disordergrid[i*self.rows+j][0]
                            jr = self.disordergrid[i*self.rows+j][1]
                            if self.grid[x,y] != 0 and (np.sqrt((x+xr-i-ir)**2+(y+yr-j-jr)**2) <= radiuslist[3] and np.sqrt((x+xr-i-ir)**2+(y+yr-j-jr)**2) > 0):
                                neighlist.append([x,y])
                    self.neighgrid.append(neighlist)

        if self.shape == "Random3D" or self.shape == "Random3D0":
            self.neighgrid = []
            for i in range(self.rows):
                for j in range(self.cols):
                    for k in range(self.rows):
                        neighlist = []
                        radiuslist = []
                        for x in range(self.rows):
                            for y in range(self.cols):
                                for z in range(self.rows):
                                    if self.grid[x,y,z] != 0 and (i != x or j != y or k != z):
                                        xr = self.disordergrid[x*self.rows*self.cols+y*self.rows+z][0]
                                        yr = self.disordergrid[x*self.rows*self.cols+y*self.rows+z][1]
                                        zr = self.disordergrid[x*self.rows*self.cols+y*self.rows+z][2]
                                        ir = self.disordergrid[i*self.rows*self.cols+j*self.rows+k][0]
                                        jr = self.disordergrid[i*self.rows*self.cols+j*self.rows+k][1]
                                        kr = self.disordergrid[i*self.rows*self.cols+j*self.rows+k][2]
                                        radiuslist.append(np.sqrt((x+xr-i-ir)**2+(y+yr-j-jr)**2+(z+zr-k-kr)**2))
                        radiuslist.sort()
                        for x in range(self.rows):
                            for y in range(self.cols):
                                for z in range(self.rows):
                                    xr = self.disordergrid[x*self.rows*self.cols+y*self.rows+z][0]
                                    yr = self.disordergrid[x*self.rows*self.cols+y*self.rows+z][1]
                                    zr = self.disordergrid[x*self.rows*self.cols+y*self.rows+z][2]
                                    ir = self.disordergrid[i*self.rows*self.cols+j*self.rows+k][0]
                                    jr = self.disordergrid[i*self.rows*self.cols+j*self.rows+k][1]
                                    kr = self.disordergrid[i*self.rows*self.cols+j*self.rows+k][2]
                                    if self.grid[x,y,z] != 0 and (np.sqrt((x+xr-i-ir)**2+(y+yr-j-jr)**2+(z+zr-k-kr)**2) <= radiuslist[5] and np.sqrt((x+xr-i-ir)**2+(y+yr-j-jr)**2+(z+zr-k-kr)**2) > 0):
                                        neighlist.append([x,y,z])
                        self.neighgrid.append(neighlist)
        
        self.magnetisation = float(np.sum(self.grid)) #first magnetisation value times size of lattice!
        self.magsq = self.magnetisation**2

        self.tot_e = self.total_energy()
        
        self.tot_e_sq = self.tot_e**2

        self.tot_e_list = []
        self.mag_list = []
        self.tot_e_sq_list = []
        self.magsq_list= []
    

    def print_grid(self):
        """Prints the grid"""
        print(self.grid)
    
    def second_neighbours(self, i, j):
        """Returns the second nearest neighbours"""
        leftup = (i-1, j-1)
        leftdown = ((i+1)%self.rows, j-1)
        rightup = (i-1, (j+1)%self.cols)
        rightdown = ((i+1)%self.rows, (j+1)%self.cols)
        return [leftup, leftdown, rightup, rightdown]

    def neighbours(self, i, j, k=0):
        """Finds the neighbours of point (row = i, col=j) on the (row, col) grid
            For Periodic Boundary Conditions, returns the indices of the neighbours
            For Non-Periodic Boundary Conditions, returns the value of self.grid of the neighbours, with 0 off the edges"""
        if self.n == "nn" and self.bc == True and (self.shape == "Square" or self.shape == "NiO"):
            left = (i, j-1)
            right = (i, (j+1)%self.cols)
            up = (i-1,j)
            down = ((i+1)%self.rows, j)

            return [left, right, up, down]
        
        if self.n == "fn" and self.bc == True:
            return self.neighgrid[i*self.rows+j]
        
        if self.shape == "Tri" and self.n == "nn" and self.bc == True:
            if i%2 == 0:
                left = (i, j-1)
                right = (i, (j+1)%self.cols)
                upleft = (i-1,j-1)
                upright = (i-1,j)
                downleft = ((i+1)%self.rows, j-1)
                downright = ((i+1)%self.rows, j)

            else:
                left = (i, j-1)
                right = (i, (j+1)%self.cols)
                upleft =(i-1,j)
                upright = (i-1,(j+1)%self.cols)
                downleft = ((i+1)%self.rows, j)
                downright = ((i+1)%self.rows, (j+1)%self.cols)

            return [left, right, upleft, upright, downleft, downright]

        if self.shape == "3D" and self.n == "nn" and self.bc == True:
            left = (i, j-1, k)
            right = (i, (j+1)%self.cols, k)
            up = (i-1,j, k)
            down = ((i+1)%self.rows, j, k)
            under = (i,j,k-1)
            over = (i,j,(k+1)%self.rows)
            return [left, right, up, down, under, over]

        if self.shape == "Random" or self.shape == "Random0":  #Returns 4 nearest neighbours
            return self.neighgrid[i*self.rows+j]
        
        if self.shape == "Random3D" or self.shape == "Random3D0":
            return self.neighgrid[i*self.rows*self.cols+j*self.rows+k]

        if self.n == "nn" and self.bc == False and self.shape == "Square":
            left = (i, j-1)
            right = (i, (j+1))
            up = (i-1,j)
            down = ((i+1), j)
            
            if i == 0 and j != 0 and j != self.cols - 1: #top row
                return [self.grid[left[0], left[1]],
                        self.grid[right[0], right[1]],
                        0,
                        self.grid[down[0], down[1]]]
            
            if i == self.rows - 1 and j != 0 and j != self.cols - 1: #bottom row
                return [self.grid[left[0], left[1]],
                        self.grid[right[0], right[1]],
                        self.grid[up[0], up[1]],
                        0]
            
            if j == 0 and i != 0 and i != self.rows - 1: #left column
                return [0,
                        self.grid[right[0], right[1]],
                        self.grid[up[0], up[1]],
                        self.grid[down[0], down[1]]]
            
            if j == self.cols - 1 and i != 0 and i != self.rows - 1: #right column
                return [self.grid[left[0], left[1]],
                    0,
                    self.grid[up[0], up[1]],
                    self.grid[down[0], down[1]]]
            
            if i == 0 and j == self.cols - 1: #top right
                return [self.grid[left[0], left[1]],
                        0,
                        0,
                        self.grid[down[0], down[1]]]

            if i == 0 and j == 0: #top left
                return [0,
                        self.grid[right[0], right[1]],
                        0,
                        self.grid[down[0], down[1]]]

            if i == self.rows - 1 and j == 0: #bottom left
                return [0,
                        self.grid[right[0], right[1]],
                        self.grid[up[0], up[1]],
                        0]

            if i == self.rows - 1 and j == self.cols - 1: #bottom right
                return [self.grid[left[0], left[1]],
                        0,
                        self.grid[up[0], up[1]],
                        0]

            else:
                return [self.grid[left[0], left[1]],
                        self.grid[right[0], right[1]],
                        self.grid[up[0], up[1]],
                        self.grid[down[0], down[1]]]
                    
    def neighbours_notpbc(self,i,j):
        """Small function to return the actual indices of the non-periodic neighbours"""
        left = (i, j-1)
        right = (i, (j+1))
        up = (i-1,j)
        down = ((i+1), j)
        
        if i == 0 and j != 0 and j != self.cols - 1: #top row
            return [left,
                    right,
                    down]
        
        if i == self.rows - 1 and j != 0 and j != self.cols - 1: #bottom row
            return [left,
                    right,
                    up]
        
        if j == 0 and i != 0 and i != self.rows - 1: #left column
            return [right,
                    up,
                    down]
        
        if j == self.cols - 1 and i != 0 and i != self.rows - 1: #right column
            return [left,
                    up,
                    down]
        
        if i == 0 and j == self.cols - 1: #top right
            return [left,
                    down]
        
        if i == 0 and j == 0: #top left
            return [right,
                    down]
        
        if i == self.rows - 1 and j == 0: #bottom left
            return [right,
                    up]
        
        if i == self.rows - 1 and j == self.cols - 1: #bottom right
            return [left,
                    up]
        
        else:
            return [self.grid[left[0], left[1]],
                    self.grid[right[0], right[1]],
                    self.grid[up[0], up[1]],
                    self.grid[down[0], down[1]]]
                            
    def energy_change(self, i, j, k=0):
        """Calculates the change in enery if a spin is flipped
            Discrepancy in B.C. is hardcorded here"""
        if self.bc == False and self.shape == "Square":
            return 2 * self.grid[i,j] * (self.neighbours(i,j)[0]*self.jgrid[i, j, i, j-1]+self.neighbours(i,j)[1]*self.jgrid[i, j, i, (j+1)%self.cols]+self.neighbours(i,j)[2]*self.jgrid[i, j, i-1,j]+self.neighbours(i,j)[3]*self.jgrid[i, j, (i+1)%self.rows, j]) + 2 * self.h * self.grid[i,j]
        elif (self.shape == "3D" or self.shape == "Random3D" or self.shape == "Random3D0") and self.bc == True:
            tempchange = 0
            for elem in self.neighbours(i,j,k):
                tempchange += self.grid[elem[0], elem[1], elem[2]] * self.jgrid[i, j, k,elem[0], elem[1], elem[2]]
            return (2 * self.grid[i,j,k] * tempchange) + (2 * self.h * self.grid[i,j,k])
        elif self.shape == "NiO":
            tempchange = 0
            for elem in self.neighbours(i,j):
                tempchange += self.grid[elem[0], elem[1]] * 2.3
            for elem in self.second_neighbours(i,j):
                tempchange += self.grid[elem[0], elem[1]] * -21
            return (2 * self.grid[i,j] * tempchange)
        else:
            tempchange = 0
            for elem in self.neighbours(i,j):
                tempchange += self.grid[elem[0], elem[1]] * self.jgrid[i, j, elem[0], elem[1]]
            return (2 * self.grid[i,j] * tempchange) + (2 * self.h * self.grid[i,j])
    
    def total_energy(self):
        """Calculates the total energy of the ferromagnet, based on the current value of self.grid
            Discrepancy in B.C. is hardcorded here"""
        if self.bc == False and self.shape == "Square":
            tot_e = 0
            for i in range(self.rows):
                for j in range(self.cols):
                    tot_e = tot_e - self.grid[i,j] * (self.neighbours(i,j)[0]*self.jgrid[i, j, i, j-1]+self.neighbours(i,j)[1]*self.jgrid[i, j, i, (j+1)%self.cols]+self.neighbours(i,j)[2]*self.jgrid[i, j, i-1,j]+self.neighbours(i,j)[3]*self.jgrid[i, j, (i+1)%self.rows, j])
            return tot_e/2.0 - self.h * np.sum(self.grid)
        elif (self.shape == "3D" or self.shape == "Random3D" or self.shape == "Random3D0") and self.bc == True:
            tot_e = 0
            for i in range(self.rows):
                for j in range(self.cols):
                    for k in range(self.rows):
                        tempchange = 0
                        for elem in self.neighbours(i, j, k):
                            tempchange += self.grid[elem[0], elem[1],elem[2]] * self.jgrid[i, j, k,elem[0], elem[1], elem[2]]
                        tot_e = tot_e - self.grid[i,j,k] * tempchange
            return tot_e/2.0 - self.h * np.sum(self.grid)

        elif self.shape == "NiO":
            tot_e = 0
            for i in range(self.rows):
                for j in range(self.cols):
                    tempchange = 0
                    for elem in self.neighbours(i,j):
                        tempchange += self.grid[elem[0], elem[1]] * 2.3
                    for elem in self.second_neighbours(i,j):
                        tempchange += self.grid[elem[0], elem[1]] * -21
                        tot_e = tot_e - self.grid[i,j] * tempchange
            return tot_e/2.0 - self.h * np.sum(self.grid)
            
        else:
            tot_e = 0
            for i in range(self.rows):
                for j in range(self.cols):
                    tempchange = 0
                    for elem in self.neighbours(i, j):
                        tempchange += self.grid[elem[0], elem[1]] * self.jgrid[i, j, elem[0], elem[1]]
                    tot_e = tot_e - self.grid[i,j] * tempchange
            return tot_e/2.0 - self.h * np.sum(self.grid)

    def metropolis(self):
        """Applies the Metropolis algorithm to calculate the next grid (Simplified)"""
        if self.shape == "Square" or self.shape == "Tri" or self.shape == "NiO":
            for i in range(self.rows):
                for j in range(self.cols):
                
                    e = self.energy_change(i,j)

                    if e <= 0:
                        self.grid[i,j] *= -1
                        self.magnetisation += 2 * self.grid[i,j]
                        self.magsq = self.magnetisation**2
                        if self.n == "fn":
                            self.tot_e = self.total_energy()
                        else:
                            self.tot_e += e
                        self.tot_e_sq = self.tot_e**2
                
                    elif np.exp((-1.0 * e)/self.temp) > random.random():
                        self.grid[i,j] *= -1
                        self.magnetisation += 2 * self.grid[i,j]
                        self.magsq = self.magnetisation**2
                        if self.n == "fn":
                            self.tot_e = self.total_energy()
                        else:
                            self.tot_e += e
                        self.tot_e_sq = self.tot_e**2

        if self.shape == "Random" or self.shape == "Random0":
            for i in range(self.rows):
                for j in range(self.cols):
                    if self.grid[i,j]!=0:
                        
                        e = self.energy_change(i,j)
                        
                        if e <= 0:
                            self.grid[i,j] *= -1
                            self.magnetisation += 2 * self.grid[i,j]
                            self.magsq = self.magnetisation**2
                            if self.n == "fn":
                                self.tot_e = self.total_energy()
                            else:
                                self.tot_e += e
                            self.tot_e_sq = self.tot_e**2
                                
                        elif np.exp((-1.0 * e)/self.temp) > random.random():
                            self.grid[i,j] *= -1
                            self.magnetisation += 2 * self.grid[i,j]
                            self.magsq = self.magnetisation**2
                            if self.n == "fn":
                                self.tot_e = self.total_energy()
                            else:
                                self.tot_e += e
                            self.tot_e_sq = self.tot_e**2

        elif self.shape == "3D" or self.shape == "Random3D" or self.shape == "Random3D0":
            for i in range(self.rows):
                for j in range(self.cols):
                    for k in range(self.rows):
                        e = self.energy_change(i,j,k)
                        
                        if e <= 0:
                            self.grid[i,j,k] *= -1
                            self.magnetisation += 2 * self.grid[i,j,k]
                            self.magsq = self.magnetisation**2
                            if self.n == "fn":
                                self.tot_e = self.total_energy()
                            else:
                                self.tot_e += e
                            self.tot_e_sq = self.tot_e**2
                        
                        elif np.exp((-1.0 * e)/self.temp) > random.random():
                            self.grid[i,j,k] *= -1
                            self.magnetisation += 2 * self.grid[i,j,k]
                            self.magsq = self.magnetisation**2
                            if self.n == "fn":
                                self.tot_e = self.total_energy()
                            else:
                                self.tot_e += e
                            self.tot_e_sq = self.tot_e**2
    
        self.counter += 1

        self.tot_e_list.append(self.tot_e)
        self.tot_e_sq_list.append(self.tot_e_sq)
        self.mag_list.append(self.magnetisation)
        self.magsq_list.append(self.magsq)
        
        self.counter_list.append(self.counter)

    
    def metropolis_for_anim(self, *args):
        """Calls the metropolis algorithm and returns a mesh suitable for animation"""
        self.metropolis()
        self.flat_grid = np.reshape(self.grid,(1,np.product(self.grid.shape)))
        colour_array = self.Colour_Array(0)
        return colour_array
    
    def Colour_Array(self, i):
        """Defines the colours of the points to be plotted"""
        if self.shape == "3D" or self.shape == "Random3D" or self.shape == "Random3D0":
            ax = plt.axes(projection = '3d')
            colour_array = np.where(self.flat_grid[i] == 1, '#00ced1', '#fa8072')
            return ax.scatter3D(self.position[:,0], self.position[:,1], self.position[:,2], c = colour_array.T.flatten())
        
        if self.shape == "Square" or self.shape == "Tri":
            colour_array = np.where(self.flat_grid[i] == 1, '#00ced1', '#fa8072')
            return plt.scatter(self.position[:,0], self.position[:,1], c = colour_array.T.flatten())

        if self.shape == "Random" or self.shape == "Random0":
            colour_array1 = np.where(self.flat_grid[i] == 1, '#00ced1', '#fa8072')
            colour_array2 = np.where(self.flat_grid[i] == 0, '#ffffff', '0')
            colour_array = np.empty_like(colour_array1)
            for z, elem in enumerate(colour_array2):
                if colour_array2[z] == '0':
                    colour_array[z] = colour_array1[z]
                else:
                    colour_array[z] = colour_array2[z]
            return plt.scatter(self.position[:,0], self.position[:,1], c = colour_array.T.flatten())
                
        if self.shape == "NiO":
            colour_array1 = np.where(self.flat_grid[i] == 1, '#00ced1', '#fa8072')
            colour_array2 = np.where(self.flat_grid[i] == 0, '#bebebe', '0')
            colour_array = np.empty_like(colour_array1)
            for z, elem in enumerate(colour_array2):
                if colour_array2[z] == '0':
                    colour_array[z] = colour_array1[z]
                else:
                    colour_array[z] = colour_array2[z]
            return plt.scatter(self.position[:,0], self.position[:,1], c = colour_array.T.flatten())

    def show_neighbours(self, i, j, k=0):
        """Function to plot a point of interest and it's neighbours, based on the definition of the lattice"""
        l=0
        if (self.shape == "Square" or self.shape == "Tri"):
            if self.bc == True:
                neighlist = self.neighbours(i,j)
            if self.bc == False:
                neighlist = self.neighbours_notpbc(i,j)
            self.flat_grid = np.reshape(self.grid,(1,np.product(self.grid.shape)))
            fig = plt.figure(figsize=(8,8),dpi = 80)
            colour_array = np.where(self.flat_grid[l] == 1, '#00ced1', '#fa8072')
            plt.scatter(self.position[:,0], self.position[:,1], c = colour_array.T.flatten())
            for elem in neighlist:
                if elem != 0:
                    plt.scatter(self.position[elem[0]*self.rows+elem[1],0], self.position[elem[0]*self.rows+elem[1],1], c = '#ffff00', s=50)
            plt.scatter(self.position[i*self.rows+j,0],self.position[i*self.rows+j,1], c = '#ffff00', s=50)
            plt.title("Neighbours of the Point (%i, %i) \n $\sigma_j = 1 \longrightarrow Blue, \sigma_j = -1 \longrightarrow Red$" %(i,j))
            plt.show()
    
        if self.shape == "Random" or self.shape == "Random0":
            neighlist = self.neighbours(i,j)
            self.flat_grid = np.reshape(self.grid,(1,np.product(self.grid.shape)))
            fig = plt.figure(figsize=(8,8),dpi = 80)
            colour_array1 = np.where(self.flat_grid[l] == 1, '#00ced1', '#fa8072')
            colour_array2 = np.where(self.flat_grid[l] == 0, '#ffffff', '0')
            colour_array = np.empty_like(colour_array1)
            for z, elem in enumerate(colour_array2):
                if colour_array2[z] == '0':
                    colour_array[z] = colour_array1[z]
                else:
                    colour_array[z] = colour_array2[z]
            plt.scatter(self.position[:,0], self.position[:,1], c = colour_array.T.flatten())
            for elem in neighlist:
                plt.scatter(self.position[elem[0]*self.rows+elem[1],0], self.position[elem[0]*self.rows+elem[1],1], c = '#ffff00', s=50)
            plt.scatter(self.position[i*self.rows+j,0],self.position[i*self.rows+j,1], c = '#000000', s=50)
            plt.title("Neighbours of the Point (%i, %i) \n $\sigma_j = 1 \longrightarrow Blue, \sigma_j = -1 \longrightarrow Red$" %(i,j))
            plt.show()

        if self.shape == "NiO":
            neighlist = self.neighbours(i,j)
            self.flat_grid = np.reshape(self.grid,(1,np.product(self.grid.shape)))
            fig = plt.figure(figsize=(8,8),dpi = 80)
            colour_array1 = np.where(self.flat_grid[l] == 1, '#00ced1', '#fa8072')
            colour_array2 = np.where(self.flat_grid[l] == 0, '#bebebe', '0')
            colour_array = np.empty_like(colour_array1)
            for z, elem in enumerate(colour_array2):
                if colour_array2[z] == '0':
                    colour_array[z] = colour_array1[z]
                else:
                    colour_array[z] = colour_array2[z]
            plt.scatter(self.position[:,0], self.position[:,1], c = colour_array.T.flatten())
            for elem in neighlist:
                plt.scatter(self.position[elem[0]*self.rows+elem[1],0], self.position[elem[0]*self.rows+elem[1],1], c = '#ffff00', s=50)
            for elem in self.second_neighbours(i,j):
                plt.scatter(self.position[elem[0]*self.rows+elem[1],0], self.position[elem[0]*self.rows+elem[1],1], c = '#ffff00', s=50)
            plt.scatter(self.position[i*self.rows+j,0],self.position[i*self.rows+j,1], c = '#000000', s=50)
            plt.title("Neighbours of the Point (%i, %i) \n $\sigma_j = 1 \longrightarrow Blue, \sigma_j = -1 \longrightarrow Red$, O atoms $\longrightarrow$ Grey" %(i,j))
            plt.show()

        if self.shape == "Random3D" or self.shape == "Random3D0":
            neighlist = self.neighbours(i,j,k)
            self.flat_grid = np.reshape(self.grid,(1,np.product(self.grid.shape)))
            fig = plt.figure(figsize=(8,8),dpi = 80)
            ax = plt.axes(projection = '3d')
            colour_array1 = np.where(self.flat_grid[l] == 1, '#00ced1', '#fa8072')
            colour_array2 = np.where(self.flat_grid[l] == 0, '#ffffff', '0')
            colour_array = np.empty_like(colour_array1)
            for z, elem in enumerate(colour_array2):
                if colour_array2[z] == '0':
                    colour_array[z] = colour_array1[z]
                else:
                    colour_array[z] = colour_array2[z]
            ax.scatter3D(self.position[:,0], self.position[:,1], self.position[:,2], c = colour_array.T.flatten())
            for elem in neighlist:
                ax.scatter3D(self.position[elem[0]*self.rows*self.cols+elem[1]*self.rows+elem[2],0], self.position[elem[0]*self.rows*self.cols+elem[1]*self.rows+elem[2],1], self.position[elem[0]*self.rows*self.cols+elem[1]*self.rows+elem[2],2], c = '#000000', s=50)
            ax.scatter3D(self.position[i*self.rows*self.cols+j*self.rows+k,0],self.position[i*self.rows*self.cols+j*self.rows+k,1], self.position[i*self.rows*self.cols+j*self.rows+k,2], c = '#ffff00', s=50)
            plt.title("Neighbours of the Point (%i, %i, %i); $\sigma_j = 1 \longrightarrow Blue, \sigma_j = -1 \longrightarrow Red$" %(i,j,k))
            plt.show()

        if self.shape == "3D":
            neighlist = self.neighbours(i,j,k)
            self.flat_grid = np.reshape(self.grid,(1,np.product(self.grid.shape)))
            fig = plt.figure(figsize=(8,8),dpi = 80)
            ax = plt.axes(projection = '3d')
            colour_array = np.where(self.flat_grid[l] == 1, '#00ced1', '#fa8072')
            ax.scatter3D(self.position[:,0], self.position[:,1], self.position[:,2], c = colour_array.T.flatten())
            for elem in neighlist:
                ax.scatter3D(self.position[elem[0]*self.rows*self.cols+elem[1]*self.rows+elem[2],0], self.position[elem[0]*self.rows*self.cols+elem[1]*self.rows+elem[2],1], self.position[elem[0]*self.rows*self.cols+elem[1]*self.rows+elem[2],2], c = '#000000', s=50)
            ax.scatter3D(self.position[i*self.rows*self.cols+j*self.rows+k,0],self.position[i*self.rows*self.cols+j*self.rows+k,1], self.position[i*self.rows*self.cols+j*self.rows+k,2], c = '#fa8072', s=50)
            plt.title("Neighbours of the Point (%i, %i, %i); $\sigma_j = 1 \longrightarrow Blue, \sigma_j = -1 \longrightarrow Red$" %(i,j,k))
            plt.show()


