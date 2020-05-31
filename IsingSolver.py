#!/usr/bin/python
#import matplotlib
#matplotlib.use('TKAgg')
import matplotlib.pylab as plt
import matplotlib.animation as animation
import numpy as np
#from scipy.constants import k
#import math
import random
from scipy.ndimage import gaussian_filter

class IsingSolver:
    """A class for analysing the Ising Model"""

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
                    if (i*self.rows+j)%2 == 0:
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
            if self.shape == "3D":
                self.jgrid = -1*np.ones([self.rows, self.cols, self.rows, self.cols, self.cols, self.rows])
            else:
                self.jgrid = -1*np.ones([self.rows, self.cols, self.rows, self.cols])

        if self.jtype == 5: #Default set-up for jgrid
            if self.shape == "3D":
                self.jgrid = 5*np.ones([self.rows, self.cols, self.rows, self.rows, self.cols, self.rows])
            else:
                self.jgrid = 5*np.ones([self.rows, self.cols, self.rows, self.cols])

        if self.jtype == -5:
            if self.shape == "3D":
                self.jgrid = -5*np.ones([self.rows, self.cols, self.rows, self.cols, self.cols, self.rows])
            else:
                self.jgrid = -5*np.ones([self.rows, self.cols, self.rows, self.cols])

        if self.jtype == 10: #Default set-up for jgrid
            if self.shape == "3D":
                self.jgrid = 10*np.ones([self.rows, self.cols, self.rows, self.rows, self.cols, self.rows])
            else:
                self.jgrid = 10*np.ones([self.rows, self.cols, self.rows, self.cols])

        if self.jtype == -10:
            if self.shape == "3D":
                self.jgrid = -10*np.ones([self.rows, self.cols, self.rows, self.cols, self.cols, self.rows])
            else:
                self.jgrid = -10*np.ones([self.rows, self.cols, self.rows, self.cols])
        
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
                            self.jgrid[i, j, x, y] = random.random()*random.choice([-1,1])
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
        print self.grid
    
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
        
        if self.shape == "Tri" and self.n == "nn" and self.bc == True:    #Odd rows are shifted to the right . . . . .
            if i%2 == 0:                                                                                #     . . . . .
                left = (i, j-1)                                                                         #    . . . . .
                right = (i, (j+1)%self.cols)                                                            #     . . . . .
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
        
        if self.shape == "Random3D" or self.shape == "Random3D0": #Returns 6 nearest neighbours
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
                            
    def energy_change(self, i, j, k=0):
        """Calculates the change in enery if a spin is flipped
            Discrepancy in B.C. is hardcorded here"""
        if self.bc == False and self.shape == "Square":
            return 2 * self.grid[i,j] * (self.neighbours(i,j)[0]*self.jgrid[i, j, i, j-1]+self.neighbours(i,j)[1]*self.jgrid[i, j, i, (j+1)%self.cols]+self.neighbours(i,j)[2]*self.jgrid[i, j, i-1,j]+self.neighbours(i,j)[3]*self.jgrid[i, j, (i+1)%self.rows, j]) + 2 * self.h * self.grid[i,j]
        elif (self.shape == "3D") and self.bc == True:
            tempchange = 0
            for elem in self.neighbours(i,j, k):
                tempchange += self.grid[elem[0], elem[1], elem[2]] * self.jgrid[i, j, k,elem[0], elem[1], elem[2]]
            return (2 * self.grid[i,j,k] * tempchange) + (2 * self.h * self.grid[i,j,k])
        if self.shape == "Random0" or self.shape == "Random" or self.n == "fn":
            tempchange = 0
            for elem in self.neighbours(i,j):
                tempchange += self.grid[elem[0], elem[1]] * self.jgrid[i, j, elem[0], elem[1]]
            return (2 * self.grid[i,j] * tempchange) + (2 * self.h * self.grid[i,j])
        if self.shape == "Random3D0" or self.shape == "Random3D":
            tempchange = 0
            for elem in self.neighbours(i,j,k):
                tempchange += self.grid[elem[0], elem[1], elem[2]] * self.jgrid[i, j, k,elem[0], elem[1], elem[2]]
            return (2 * self.grid[i,j,k] * tempchange) + (2 * self.h * self.grid[i,j,k])
        if self.shape == "NiO":
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
        """Applies the Metropolis algorithm to calculate the next grid"""
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
                            self.tot_e = self.total_energy()
                            self.tot_e_sq = self.tot_e**2
                        
                        elif np.exp((-1.0 * e)/self.temp) > random.random():
                            self.grid[i,j] *= -1
                            self.magnetisation += 2 * self.grid[i,j]
                            self.magsq = self.magnetisation**2
                            self.tot_e = self.total_energy()
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
                            if self.n == "fn" or self.shape == "Random3D" or self.shape == "Random3D0":
                                self.tot_e = self.total_energy()
                            else:
                                self.tot_e += e
                            self.tot_e_sq = self.tot_e**2
                        
                        elif np.exp((-1.0 * e)/self.temp) > random.random():
                            self.grid[i,j,k] *= -1
                            self.magnetisation += 2 * self.grid[i,j,k]
                            self.magsq = self.magnetisation**2
                            if self.n == "fn" or self.shape == "Random3D" or self.shape == "Random3D0":
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
        #print self.counter
    
    def energymag_temp_start(self, warmup, N):
        """Finds the 4 different measures against temperature. Doesn't plot them, to be fed into plotting functions"""
        z = 20
        self.templist = np.linspace(0.1, 5, z)
        self.av_energy_list = np.zeros(z)
        self.av_energy_errors = np.zeros(z)
        self.av_energysq_list = np.zeros(z)
        self.av_energysq_errors = np.zeros(z)
        self.av_mag_list = np.zeros(z)
        self.av_mag_errors = np.zeros(z)
        self.av_magsq_list = np.zeros(z)
        self.av_magsq_errors = np.zeros(z)
        for n, temp in enumerate(self.templist):
            print n
            self.temp = temp
            if self.shape == "Square" or self.shape == "Tri":
                self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols))
            if self.shape == "Random":
                self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols))
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
            if self.shape == "Random0":
                self.grid = np.random.choice((-1, 0, 1), size=(self.rows, self.cols))
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
                self.dim = np.sum(np.abs(self.grid))
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
            if self.shape == "3D":
                self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols, self.rows))
            if self.shape == "NiO":
                self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols))
                for i in range(self.rows):
                    for j in range(self.cols):
                        if (i*self.rows+j)%2 == 0:
                            self.grid[i,j] = 0
            if self.n == "fn" and self.bc == True and self.shape == "Square":     #For Random Neighbours choose different random neighbours each time
                list = [1,2,3]
                self.neighgrid = []
                for i in range(self.rows):
                    for j in range(self.cols):
                        left = (i, j-random.choice(list))
                        right = (i, (j+random.choice(list))%self.cols)
                        up = (i-random.choice(list),j)
                        down = ((i+random.choice(list))%self.rows, j)
                        self.neighgrid.append([left, right, up, down])
            self.tot_e = self.total_energy()
            self.tot_e_sq = self.tot_e**2
            self.magnetisation = float(np.sum(self.grid))
            self.magsq = self.magnetisation**2
            self.counter = 0
            self.tot_e_list = []
            self.tot_e_sq_list = []
            self.mag_list = []
            self.magsq_list = []

            for j in range(warmup+N):
                self.metropolis()
            self.av_energy_list[n] = np.mean(self.tot_e_list[-N:])/(self.dim)
            self.av_energy_errors[n] = np.std(self.tot_e_list[-N:])/(self.dim)
            self.av_energysq_list[n] = np.mean(self.tot_e_sq_list[-N:])/(self.dim)**2
            self.av_energysq_errors[n] = np.std(self.tot_e_sq_list[-N:])/(self.dim)**2
            self.av_mag_list[n] = np.mean(np.abs(self.mag_list[-N:]))/(self.dim)
            self.av_mag_errors[n] = np.std(np.abs(self.mag_list[-N:]))/(self.dim)
            self.av_magsq_list[n] = np.mean(np.abs(self.magsq_list[-N:]))/(self.dim)**2
            self.av_magsq_errors[n] = np.std(np.abs(self.magsq_list[-N:]))/(self.dim)**2
    
        self.heatcap_list = (self.av_energysq_list - self.av_energy_list**2)/self.temp**2
        self.magsus_list = (self.av_magsq_list - self.av_mag_list**2)*self.temp
        self.heatcap_errors = np.sqrt(self.av_energysq_errors**2 + (2*self.av_energy_list**2*np.divide(self.av_energy_errors,self.av_energy_list))**2)/self.temp**2/100
        self.magsus_errors = np.sqrt(self.av_magsq_errors**2 + (2*self.av_mag_list**2*np.divide(self.av_mag_errors,self.av_mag_list))**2)*self.temp/100
        
        f = open("4PlotsDataRanCoupl.txt", "w+")
        f.write("Temperature List = ")
        for num in range(z):
            f.write("%1.9f, " %(self.templist[num]))
        f.write("\n")
        f.write("Av Energy List = ")
        for num in range(z):
            f.write("%1.9f, " %(self.av_energy_list[num]))
        f.write("\n")
        f.write("Magnetisation List = ")
        for num in range(z):
            f.write("%1.9f, " %(self.av_mag_list[num]))
        f.write("\n")
        f.write("Heat Capacity List = ")
        for num in range(z):
            f.write("%1.9f, " %(self.heatcap_list[num]))
        f.write("\n")
        f.write("Mag Susc List = ")
        for num in range(z):
            f.write("%1.9f, " %(self.magsus_list[num]))
        f.write("\n")
        f.write("Temperature Errors = ")
        for num in range(z):
            f.write("%1.9f, " %(self.templist[num]))
        f.write("\n")
        f.write("Av Energy Errors = ")
        for num in range(z):
            f.write("%1.9f, " %(self.av_energy_errors[num]))
        f.write("\n")
        f.write("Magnetisation Errors = ")
        for num in range(z):
            f.write("%1.9f, " %(self.av_mag_errors[num]))
        f.write("\n")
        f.write("Heat Capacity Errors = ")
        for num in range(z):
            f.write("%1.9f, " %(self.heatcap_errors[num]))
        f.write("\n")
        f.write("Mag Susc Errors = ")
        for num in range(z):
            f.write("%1.9f, " %(self.magsus_errors[num]))
        f.write("\n")
        f.close()

    def find_curie_temp(self):
        """Function to determine the Curie temperature from the array found in the start function"""
        smooth = gaussian_filter(self.av_mag_list,3.)
        dx = np.diff(self.templist)
        dy = np.diff(smooth)
        slope = np.abs(dy/dx)
        max_slope = np.abs(max(slope))
        for i,s in enumerate(slope):
            if s == max_slope:
                curie_temp = self.templist[i]
                print curie_temp
                curie_temp_index = i
        return [curie_temp, curie_temp_index]
            
    def energymag_temp_plot(self, warmup, N):
        """Module to plot the 4 measures calculated in the start function"""
        self.energymag_temp_start(warmup, N)
        curie_temp = self.find_curie_temp()[0]
        
        plt.errorbar(self.templist, self.av_energy_list, yerr=self.av_energy_errors, fmt = 'o')
        plt.axvline(curie_temp, linestyle = '--', color = 'g', label = "Tc = %1.1f J/Kb" %curie_temp)
        plt.title("Average Energy Against Temperature \n 10x10 Lattice with Random Neighbours, J = %i" %self.jtype)
        plt.xlabel("Temp [J/Kb]")
        plt.ylabel("Av. Energy")
        plt.legend(loc=2)
        plt.show()
        
        plt.errorbar(self.templist, self.av_mag_list, yerr=self.av_mag_errors, fmt = 'o')
        plt.axvline(curie_temp, linestyle = '--', color = 'g', label = "Tc = %1.1f J/Kb" %curie_temp)
        plt.title("Average Magnetisation Against Temperature \n 10x10 Lattice with Random Neighbours, J = %i" %self.jtype)
        plt.xlabel("Temp [J/Kb]")
        plt.ylabel("Av. Magnetisation")
        plt.legend()
        plt.show()
        
        plt.errorbar(self.templist, self.heatcap_list, yerr=self.heatcap_errors, fmt = 'o')
        plt.axvline(curie_temp, linestyle = '--', color = 'g', label = "Tc = %1.1f J/Kb" %curie_temp)
        plt.title("Heat Capacities Against Temperature \n 10x10 Lattice with Random Neighbours, J = %i" %self.jtype)
        plt.xlabel("Temp [J/Kb]")
        plt.ylabel("C$_v$")
        plt.legend()
        plt.show()

        plt.errorbar(self.templist, self.magsus_list, yerr=self.magsus_errors, fmt = 'o')
        plt.axvline(curie_temp, linestyle = '--', color = 'g', label = "Tc = %1.1f J/Kb" %curie_temp)
        plt.title("Average Magnetic Susceptibility Against Temperature \n 10x10 Lattice with Random Neighbours, J = %i" %self.jtype)
        plt.xlabel("Temp [J/Kb]")
        plt.ylabel("$\chi$ [Kb]")
        plt.legend()
        plt.show()

    def h_dependence(self, hlist, warmup, N):
        """Module investigates the dependence of the model on the Magnetic Field, h"""
        h_av_energy_list = []
        h_av_energy_errors = []
        h_av_mag_list = []
        h_av_mag_errors = []
        h_heatcap_list = []
        h_magsus_list = []
        h_heatcap_errors = []
        h_magsus_errors = []
        for h in hlist:
            self.h = h
            self.energymag_temp_start(warmup, N)
            h_av_energy_list.append(self.av_energy_list)
            h_av_energy_errors.append(self.av_energy_errors)
            h_av_mag_list.append(self.av_mag_list)
            h_av_mag_errors.append(self.av_mag_errors)
            h_heatcap_list.append(self.heatcap_list)
            h_magsus_list.append(self.magsus_list)
            h_heatcap_errors.append(self.heatcap_errors)
            h_magsus_errors.append(self.magsus_errors)
                
        for i, h in enumerate(hlist):
            plt.errorbar(self.templist, h_av_energy_list[i], yerr=h_av_energy_errors[i], label = "h = %i" %h, fmt = 'o')
        plt.title("Average Energy Against Temperature \n 10x10 Lattice")
        plt.xlabel("Temp [J/Kb]")
        plt.ylabel("Av. Energy")
        plt.legend()
        plt.show()

        for i, h in enumerate(hlist):
            plt.errorbar(self.templist, h_av_mag_list[i], yerr=h_av_mag_errors[i], label = "h = %i" %h, fmt = 'o')
        plt.title("Average Magnetisation Against Temperature \n 10x10 Lattice")
        plt.xlabel("Temp [J/Kb]")
        plt.ylabel("Av. Magnetisation")
        plt.legend()
        plt.show()
            
        for i, h in enumerate(hlist):
            plt.errorbar(self.templist, h_heatcap_list[i], yerr=h_heatcap_errors[i], label = "h = %i" %h, fmt = 'o')
        plt.title("Heat Capacities Against Temperature \n 10x10 Lattice")
        plt.xlabel("Temp [J/Kb]")
        plt.ylabel("C$_v$")
        plt.legend()
        plt.show()

        for i, h in enumerate(hlist):
            plt.errorbar(self.templist, h_magsus_list[i], yerr=h_magsus_errors[i], label = "h = %i" %h, fmt = 'o')
        plt.title("Average Magnetic Susceptibility Against Temperature \n 10x10 Lattice")
        plt.xlabel("Temp [J/Kb]")
        plt.ylabel("$\chi$ [Kb]")
        plt.legend()
        plt.show()

    def hysteresis_loop(self, A, B):
        """Plots the hysteresis loop, i.e. plots the average magnetisation as we increase h and then decrease h, i.e. increase and decrease the magnetic field"""
        hlist1 = np.linspace(0.001, B, A)
        hlist2 = np.linspace(B, -B, 2*A)
        hlist3 = np.linspace(-B, B, 2*A)
        mags1 = np.zeros(A)
        mags2 = np.zeros(2*A)
        mags3 = np.zeros(2*A)
        for i, h in enumerate(hlist1):
            self.h = h
            mags1[i] = self.magnetisation/self.dim
            self.metropolis()
        for i, h in enumerate(hlist2):
            self.h = h
            mags2[i] = self.magnetisation/self.dim
            self.metropolis()
        for i, h in enumerate(hlist3):
            self.h = h
            mags3[i] = self.magnetisation/self.dim
            self.metropolis()

        smooth2 = gaussian_filter(mags2, 3.)
        smooth3 = gaussian_filter(mags3, 3.)
        dec_h = np.mean(hlist2[np.where(np.round(smooth2, 0) == 0)])
        inc_h = np.mean(hlist3[np.where(np.round(smooth3, 0) == 0)])

        plt.plot(hlist1, mags1, label="Warm Up")
        plt.plot(hlist2, mags2, label="Decreasing h")
        txt1 = (" %1.1f" %dec_h)
        plt.text(dec_h, 0, txt1)
        plt.scatter(dec_h, 0)
        plt.plot(hlist3, mags3, label="Increasing h")
        txt2 = (" %1.1f" %inc_h)
        plt.text(inc_h, 0, txt2)
        plt.scatter(inc_h, 0)
        plt.title("Hysteresis Loop for a 10x10 Square Lattice with J = 1")
        plt.ylabel("Av. Magnetisation")
        plt.xlabel("Strength of External Magnetic Field, h")
        plt.ylim(-1.1,1.1)
        plt.legend(loc=2)
        plt.grid()
        plt.show()

    def hysteresis_J_dep(self, A, B, jlist):
        """Plots hysteresis loops for J equal to the input list"""
        hlist1 = np.linspace(0, B, A)
        hlist2 = np.linspace(B, -B, 2*A)
        hlist3 = np.linspace(-B, B, 2*A)
        hlist = np.concatenate([hlist2,hlist3])
        for j in jlist:
            self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols))
            self.magnetisation = float(np.sum(self.grid))
            self.jgrid = j*np.ones([self.rows, self.cols, self.rows, self.cols])
            mags1 = np.zeros(A)
            mags2 = np.zeros(2*A)
            mags3 = np.zeros(2*A)
            for i, h in enumerate(hlist1):
                self.h = h
                mags1[i] = self.magnetisation/self.dim
                self.metropolis()
            for i, h in enumerate(hlist2):
                self.h = h
                mags2[i] = self.magnetisation/self.dim
                self.metropolis()
            for i, h in enumerate(hlist3):
                self.h = h
                mags3[i] = self.magnetisation/self.dim
                self.metropolis()
        
            smooth2 = gaussian_filter(mags2, 3.)
            smooth3 = gaussian_filter(mags3, 3.)
            dec_h = np.mean(hlist2[np.where(np.round(smooth2, 0) == 0)])
            inc_h = np.mean(hlist3[np.where(np.round(smooth3, 0) == 0)])

            mags = np.concatenate([mags2,mags3])

            plt.plot(hlist, mags, label = "J = %i" %j)

            txt1 = (" %1.1f" %dec_h)
            plt.text(dec_h, 0, txt1)
            plt.scatter(dec_h, 0)
            txt2 = (" %1.1f" %inc_h)
            plt.text(inc_h, 0, txt2)
            plt.scatter(inc_h, 0)

        plt.title("Hysteresis Loop for a 10x10 Square Lattice for Various J")
        plt.ylabel("Av. Magnetisation")
        plt.xlabel("Strength of External Magnetic Field, h")
        plt.ylim(-1.1,1.1)
        plt.legend(loc=2)
        plt.grid()
        plt.show()

    def iterations_for_convergence_plot(self):
        """Plots the average number of iterations for convergence for the different sizes of a lattice"""
        if self.shape == "Square" and self.jtype == 1:
            z = 10   #number of lattices to average over
            dimlist = [10,20,30,40,50,60,70,80,90,100]
            av_conv_list = np.zeros(10)
            av_conv_errors = np.zeros(10)
            f = open("iterations.txt", "w+")
            f.close()
            for m, dim in enumerate(dimlist):
                counter_list = np.zeros(z)
                self.rows = dim
                self.cols = dim
                self.jgrid = np.ones([self.rows, self.cols, self.rows, self.cols])
                for n in range(z):
                    print "Dimension = ", dim
                    print "Lattice Number = ", n
                    self.grid = np.random.choice((-1, 1), size=(self.rows, self.cols))
                    self.magnetisation = float(np.sum(self.grid))
                    self.counter = 0
                    while np.abs(self.magnetisation/(self.cols*self.rows)) != 1 and self.counter < 5000:
                        self.metropolis()
                    if self.counter < 5000:
                        counter_list[n] = self.counter
                    elif n == 0:
                       n = 0
                    else:
                        counter_list[n] = counter_list[n-1]
                f = open("iterations.txt", "a+")
                f.write("%i \n" %dim)
                for num in range(z):
                    f.write("%i, " %(counter_list[num]))
                f.write("\n")
                f.close()
                av_conv_list[m] = np.mean(counter_list)
                av_conv_errors[m] = np.std(counter_list)
            
            plt.scatter(dimlist, av_conv_list)
            plt.errorbar(dimlist, av_conv_list, yerr = av_conv_errors, fmt = None)
            plt.title("Average Number of Iterations to Converge Against Size \n For Square Model with T = %1.1f [J/Kb]" %self.temp)
            plt.xlabel("Length of Square Side")
            plt.ylabel("Av. # of Iterations")
            plt.show()


    def plot_energy(self, N, type):
        """Plots the energy at every iteration of Metropolis algorithm up to N steps"""
        if type == "Total":
            for i in range(N):
                #print i
                self.metropolis()
            plt.plot(self.counter_list, self.tot_e_list)
            plt.title("Total Energy of Ising Model")
            plt.ylabel("Total Energy")
            plt.xlabel("Iteration Number")
            plt.show()
        if type == "Average":
            for i in range(N):
                #print i
                self.metropolis()
            plt.plot(self.counter_list, self.av_e_list)
            plt.title("Average Energy of Ising Model")
            plt.ylabel("Average Energy")
            plt.xlabel("Iteration Number")
            plt.show()

    def plot_magnetisation(self, N):
        """Plots the magnetisation at every iteration of Metropolis algorithm up to N steps"""
        for i in range(N):
            self.metropolis()
        plt.plot(self.counter_list, self.av_mag_list)
        plt.scatter(0, 1.1, color="white")
        plt.title("Average Magnetisation of Ising Model \n %i x %i Grid, T = %1.1f J/Kb" %(self.rows, self.cols, self.temp))
        plt.ylabel("Average Magnetisation")
        plt.xlabel("Iteration Number")
        plt.show()



