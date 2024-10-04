import numpy as np
from scipy.sparse import coo_matrix
import cvxopt
import cvxopt.cholmod
from gridMesher import GridMesh

class FE:
    def __init__(self, mesh, matProp, bc):
        self.mesh = GridMesh(mesh, matProp, bc)

    def solve(self, density):
        self.u=np.zeros((self.mesh.ndof,1))
        E = self.mesh.material['E']*(1.0e-6 + density)**self.mesh.material['penal']
        sK = np.einsum('i,ijk->ijk',E, self.mesh.KE).flatten()
        K = coo_matrix((sK,(self.mesh.iK,self.mesh.jK)),shape=(self.mesh.ndof,self.mesh.ndof)).tocsc()
        K = self.deleterowcol(K,self.mesh.fixed,self.mesh.fixed).tocoo()
        K = cvxopt.spmatrix(K.data,K.row.astype(int),K.col.astype(int))
        B = cvxopt.matrix(self.mesh.f[self.mesh.free,0])
        cvxopt.cholmod.linsolve(K,B)
        self.u[self.mesh.free,0]=np.array(B)[:,0]
        uElem = self.u[self.mesh.edofMat].reshape(self.mesh.numElems,self.mesh.numDOFPerElem)
        self.Jelem  = np.einsum('ik,ik->i',np.einsum('ij,ijk->ik',uElem, self.mesh.KE),uElem)
        return self.u, self.Jelem
    def deleterowcol(self, A, delrow, delcol):
        m = A.shape[0]
        keep = np.delete(np.arange(0, m), delrow)
        A = A[keep, :]
        keep = np.delete(np.arange(0, m), delcol)
        A = A[:, keep]
        return A


