#!/usr/bin/env python
# coding: utf-8

#Code used to estimate critical exponents Beta and gamma

import numpy as np
import random
import scipy

class lattice():
    "Represents a lattice"

    def __init__(self, L):
        self.L = L
        self.N = self.L * self.L
        self.Empty = -self.N - 1
        self.Boundary = -self.N - 2
        self.ptr = np.zeros(self.N, dtype=int)
        self.nn = np.zeros((self.N, 4), dtype=int)
        self.order = np.zeros(self.N, dtype=int)

    def nearestneighbour(self):
        """constructs a numpy.ndarray of the nearest neighbours of each site in a square lattice
           with N=L*L sites."""

        for i in range(self.N):
            self.nn[i, 0] = (i + 1) % self.N
            self.nn[i, 1] = (i - 1) % self.N
            self.nn[i, 2] = (i + self.L) % self.N
            self.nn[i, 3] = (i - self.L) % self.N

            if i % self.L == 0:
                self.nn[i, 1] = self.Boundary

            if (i + 1) % self.L == 0:
                self.nn[i, 0] = self.Boundary

            if np.floor(i / self.L) == 0:
                self.nn[i, 3] = self.Boundary

            if np.floor(i / self.L) == self.L - 1:
                self.nn[i, 2] = self.Boundary

    def permutate(self):
        """constructs a random permutation of the site labels (integers from 0 to N-1) using the algorithm
        given earlier and stores it in a numpy.ndarray """

        for i in range(self.N):
            self.order[i] = i

        for i in range(self.N):
            # generates a random number j between i and N
            j = random.randint(i, self.N - 1)

            # swaps the positions of sites i and j
            temp = self.order[i]
            self.order[i] = self.order[j]
            self.order[j] = temp

    def findroot(self, i):
        """Recursive function that returns the label (an integer from 0 to N-1) of the root site of the cluster
        that the input site (label i) is part of. Also carries out path compression by
        setting the value of ptr[i] to be the (label of the) root of the tree that site i is part of . """
        step = 0
        
        if self.ptr[i] < 0:
            return i, step
        else:
            self.ptr[i] = self.findroot(self.ptr[i])[0]
            step += 1
            return self.ptr[i], step

    def percolate(self):
        """Model for site percolation on a square lattice with N sites using the Newman-Ziff
        algorithim."""

        n = []
        largest_cluster_size = []
        average_cluster_size = []
        largest = 0  # size of largest cluster
        total = 0  # total number of clusters
        steps = 0 #number of steps taken through trees
        sum_squared = np.zeros(self.N, dtype='int64')

        for i in range(self.N):
            # constructs lattice with occupied sites along two (opposite) edges and all other sites unoccupied
            self.ptr[i] = self.Empty

        for i in range(self.N):
            # occupies sites in order
            root_1 = site_1 = self.order[i]
            self.ptr[site_1] = -1
            total += 1
            sum_squared[i] += 1

            for j in range(4):
                # finds nearest neighbours
                site_2 = self.nn[site_1, j]

                if site_2 != self.Boundary:

                    if self.ptr[site_2] != self.Empty:
                        # If nearest neighbour is occupied, find root of cluster it is part of
                        root_2 = self.findroot(site_2)[0]
                        steps += self.findroot(site_2)[1]

                        if root_2 != root_1:
                            # If root of site and its nearest neighbour are different
                            # adds smaller cluster to larger cluster
                            total -= 1
                            sum_squared[i] -= (self.ptr[root_2])**2
                            sum_squared[i] -= (self.ptr[root_1])**2

                            if self.ptr[root_1] > self.ptr[root_2]:
                                self.ptr[root_2] += self.ptr[root_1]
                                self.ptr[root_1] = root_2
                                root_1 = root_2
                                sum_squared[i] += (self.ptr[root_2])**2

                            else:
                                self.ptr[root_1] += self.ptr[root_2]
                                self.ptr[root_2] = root_1
                                sum_squared[i] += (self.ptr[root_1])**2

                            if -self.ptr[root_1] > largest:
                                # Keeps track of largest cluster in lattice
                                largest = -self.ptr[root_1]

                
            if i < self.N -1:
                sum_squared[i+1] = sum_squared[i] 

            n.append(i + 1)
            largest_cluster_size.append(largest)
            average_cluster_size.append(sum_squared[i]/(i+1))

        return largest_cluster_size, average_cluster_size
    

Lattice_sizes = [10, 30, 50, 70, 90, 100, 120, 150, 200, 250, 300, 400, 500]

#Find percolation strength, P_inf, at p=p_c for different lattice sizes
#Find average cluster size, S, at p=p_c for different lattice sizes

P_inf = [] #Percolation strength
S = [] #Average cluster size

for i in Lattice_sizes:
    sample = lattice(i)
    p = []
    s = []
    for k in range(10):
        L1 = np.zeros(i ** 2)
        L2 = np.zeros(i ** 2)
        P1 = 0
        P2 = 0
        for j in range(1000):
            #Find average cluster size and largest cluster size at each step
            #Average over 1000 runs
            lattice.nearestneighbour(sample)
            lattice.permutate(sample)
            largest, average_cluster = lattice.percolate(sample)

            L1 += np.array(largest)/1000
            L2 += np.array(average_cluster)/1000


        for n in range(i**2):
            #perform convolution with Binomial distribution at p=p_c
            #to give average cluster size and largest cluster size at p=p_c
            P1 += (L1[n]/i**2)*(scipy.stats.binom.pmf(n+1, i**2, 0.592746))
            P2 += (L2[n]) * (scipy.stats.binom.pmf(n+1, i**2, 0.592746))
        
       
        s.append(P2)
        p.append(P1)
        
    #Repeat 10 times for each lattice size and average
    S.append(np.mean(s))
    P_inf.append(np.mean(p))

data2 = np.array([Lattice_sizes, P_inf])
data3 = np.array([Lattice_sizes, S])
np.save('Percolationstrength.npy', data2)
np.svae('averagecluster.npy', data3)





   
     

        
    
    
        




