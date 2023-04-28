#!/usr/bin/env python
# coding: utf-8

# In[ ]:

# Code used to test the performance of the Newmann - Ziff algorithm

import numpy as np
import random
import scipy
import timeit


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

        Percolation = np.zeros(self.N, dtype=int)
        n = []
        largest_cluster= []
        average_cluster = []
        largest = 0  # size of largest cluster
        total = 2  # total number of clusters
        steps = 0 #number of steps taken through trees
        sum_squared = np.zeros(self.N, dtype='int64')
        self.ptr[0] = -1
        self.ptr[self.L - 1] = -1

        for i in range(self.N):
            # constructs lattice with occupied sites along two (opposite) edges and all other sites unoccupied
            self.ptr[i] = self.Empty

            if i % self.L == 0 and i != 0:
                self.ptr[i] = 0
                self.ptr[0] -= 1

            if (i + 1) % self.L == 0 and i != self.L - 1:
                self.ptr[i] = self.L - 1
                self.ptr[self.L - 1] -= 1

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

            if self.findroot(0)[0] == self.findroot(self.L - 1)[0]:
                # if two initial clusters have same root, percolation has occurred
                Percolation[i] = 1
                
            if i < self.N -1:
                sum_squared[i+1] = sum_squared[i] 

            n.append(i + 1)
            largest_cluster.append(largest)
            average_cluster.append(sum_squared[i]/(i+1))

        return n, Percolation, largest_cluster, average_cluster, steps

Lattice_sizes = np.array([20, 50, 80, 100, 120, 150, 200, 250, 300])

# find the number of steps taken through trees and the time taken during each run of the algorithm for different lattice sizes

steps = []
times = []

for i in Lattice_sizes:
    sample = lattice(i)
    s = []
    t = []
    for j in range(100):
        start = timeit.default_timer()
        lattice.nearestneighbour(sample)
        lattice.permutate(sample)
        s.append(lattice.percolate(sample)[4])
        stop = timeit.default_timer()
        t.append(stop - start)
    
    times.append(np.mean(t))
    steps.append(np.mean(s))

data4 = np.array([Lattice_sizes**2, steps])
data5 = np.array([Lattice_sizes**2, times])
np.save('totalsteps.npy', data4)
np.save('times.npy', data5)


