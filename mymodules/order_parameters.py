import numpy as np
import scipy.special as sp


class Descriptor:

    def __init__(self, coord, box, type_list = None):
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
        # array with the distances of the neighbours in order and list of neighbours in order
        self.nn_dist, self.nn_list = self.compute_nn(dist)
        # array with all the relative positions of the neighbours with respect to atom iatom
        self.nn_matrixdiff = dist[self.nn_list]

    def compute_q(self, iatom, itype = None):

        self.compute_descr(iatom, itype)

        # take just the first four nearest neighbour
        q = 0
        for i in range(3):
            for j in range(i+1,4):
                dot_prod = np.dot(self.nn_matrixdiff[i+1], self.nn_matrixdiff[j+1])
                i_norm = np.linalg.norm(self.nn_matrixdiff[i+1])
                j_norm = np.linalg.norm(self.nn_matrixdiff[j+1])
                cosphi = dot_prod/i_norm/j_norm
                q += (cosphi + 1/3)**2
        q = 1 - 3/8*q
        return q

    def apply_pbc(self, dist):
        box = np.tile(self.box,(np.shape(dist)[0],1))
        dist -= np.rint(dist/box) * box
        return dist

    def compute_nn(self, dist):
        dist = np.linalg.norm(dist, axis = 1)
        ord = dist.argsort()
        return dist[ord[:]] , ord[:]
