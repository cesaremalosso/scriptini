import numpy as np
import scipy.special as sp


class descriptor:

    def __init__(self, coord, natoms, box, type_list = None):
        self.natoms = natoms
        self.coord = coord
        self.box = box
        self.type_list = type_list

    def compute_descr(self, iatom, itype = None):

        # compute distances with respect to iatom
        if not itype:
            dist = self.coord - self.coord[iatom]
            dist = self.apply_pbc(dist)
        else:
            ss_list = []
            for jatom, ntype in enumerate(self.type_list):
                if ntype == itype:
                    ss_list.append(jatom)
            dist = self.coord[ss_list] - self.coord[iatom]
            dist = self.apply_pbc(dist)
        self.nn_dist, self.nn_list = self.compute_nn(dist)

    def compute_q(self, iatom, itype = None):

        if not self.nn_dist:
            self.compute_descr(iatom, itype)

        # take just the first four nearest neighbour
        q = 0
        for i in range(3):
            for j in range(i+1):
                q += (np.dot(self.nn_dist[i], self.nn_dist[j]) + 1/3)**2
        q = 1 - 3/8*q
        return q

    def apply_pbc(self, dist):
        box = np.tile(self.box,(np.shape(dist)[0],1))
        dist -= np.rint(dist/box) * box
        return dist

    def compute_nn(self, dist):
        dist = np.linalg.norm(dist, axis = 1)
        ord = dist.argsort()
        return dist[ord[1:]] , ord[1:]
