# Versions
# Numpy 1.18.1
# Pytorch 1.5.0
# scipy 1.4.1
# cvxopt 1.2.0
import time
from FE import FE
from lib.workspace import *
from lib.models.decoder import *
from lib.utils import *
import matlab.engine
eng = matlab.engine.start_matlab()
from tqdm import tqdm
from sklearn.neighbors import NearestNeighbors

def to_np(x):
    return x.detach().cpu().numpy()
class TopologyOptimizer:
    # -----------------------------#
    def __init__(self, mesh, matProp, bc, desiredVolumeFraction, densityProjection, res, overrideGPU=True):
        self.device = self.setDevice(overrideGPU)
        self.FE = FE(mesh, matProp, bc)
        self.desiredVolumeFraction = desiredVolumeFraction
        self.density = self.desiredVolumeFraction * np.ones((self.FE.mesh.numElems))
        self.densityProjection = densityProjection
        self.resolution = res
        self.len = 1.0
        c = self.generate_points(self.resolution, self.len )
        xy = torch.tensor(c, dtype=torch.float32)
        self.xy_e = xy.clone().detach().requires_grad_(True).float().view(-1, 3).to(self.device)
        self.objective = 0.0

        II, JJ, SS = self.generate_filterM()
        II = torch.tensor(II)
        II = II.int()
        JJ = torch.tensor(JJ)
        JJ = JJ.int()
        SS = torch.tensor(SS)
        SS = torch.reshape(SS, (-1,))
        ind = torch.cat((torch.reshape(II, (1, -1)), torch.reshape(JJ, (1, -1))), dim=0)
        self.LL = torch.sparse_coo_tensor(indices=ind, values=SS, size=(self.FE.mesh.numElems, self.FE.mesh.numElems), dtype=torch.float32)

        ###########################################
        experiment_directory_source = '../Morphology Learning/experiments/cellculture'
        specs_filename = os.path.join(experiment_directory_source, "specs.json")
        specs = json.load(open(specs_filename))
        latent_size = specs["CodeLength"]
        decoder_source = DeepSDF(latent_size, **specs["NetworkSpecs"])
        decoder_source = torch.nn.DataParallel(decoder_source)
        saved_model_state = torch.load(os.path.join(experiment_directory_source, model_params_subdir, "latest-best.pth"))
        decoder_source.load_state_dict(saved_model_state["model_state_dict"])
        self.decoder = decoder_source.module.cuda()
        reconstruction_codes_dir = os.path.join(experiment_directory_source, latent_codes_subdir)
        latent_filename = os.path.join(reconstruction_codes_dir, "latest-best.pth")
        self.latent = torch.load(latent_filename)["latent_codes"]["weight"]
    # -----------------------------#
    def setDevice(self, overrideGPU):
        if (torch.cuda.is_available() and (overrideGPU == False)):
            device = torch.device("cuda:0")
            print("GPU enabled")
        else:
            device = torch.device("cpu")
            print("Running on CPU")
        return device

    def generate_points(self, N, l):
        step = l / N
        v = np.zeros((N * N * N, 3))
        cen = np.array([-l / 2 + step / 2, l / 2 - step / 2, -l / 2 + step / 2])
        k = 0
        for z in range(N):
            for x in range(N):
                for y in range(N):
                    v[k, :] = cen + np.array([x, -(y), z]) * step
                    k = k + 1
        return v

    def generate_filterM(self,):
        p = self.FE.mesh.elemCenters
        kk = 64
        nbrs = NearestNeighbors(n_neighbors=kk, algorithm='ball_tree').fit(p)
        dis, ind = nbrs.kneighbors(p)
        dis1 = np.exp(-2 * dis)
        dd = np.sum(dis1, 1)
        dis1 = dis1 / np.reshape(dd, (-1, 1))
        a = np.reshape(np.arange(0, p.shape[0]), (-1, 1))
        b = np.tile(a, (1, kk))
        II = b.flatten()
        JJ = ind.flatten()
        SS = dis1.flatten()

        return II, JJ, SS
    # -----------------------------#
    def projectDensity(self, x):
        if (self.densityProjection['isOn']):
            b = self.densityProjection['sharpness']
            xp = (torch.exp(2*x) - 1)*(torch.exp(-2*x) + 1)
            xp = ((np.tanh(0.0 * b) + torch.tanh(b * (xp - 0.0))) / (np.tanh(0.0 * b) + np.tanh(b*(1-0.0))) + 1) / 2
        return xp

    def densityFilter(self, x):
        xt = torch.matmul(self.LL, torch.reshape(x, (-1, 1)))
        xt = torch.reshape(xt, (-1, ))
        return xt

    def structural_opt(self, maxepoch):
        x = self.latent[0]
        x.requires_grad = True
        lr = 5e-3
        self.decoder.eval()
        compliance = []
        vol = []
        alphaMax = 500
        alphaIncrement = 1
        alpha = 10
        self.optimizer = torch.optim.Adam([x], lr=lr)
        print("Starting optimization:")

        for e in tqdm(range(maxepoch)):
            self.optimizer.zero_grad()
            sdf_e = torch.ones(self.FE.mesh.numElems)
            sdf_e[0:self.FE.mesh.nelx*self.FE.mesh.nely*(self.FE.mesh.nelz-1)] = \
                torch.flatten(decode_sdf(self.decoder, x, self.xy_e).to(self.device))
            rhoFilter = self.densityFilter(sdf_e)
            rhoElem = self.projectDensity(rhoFilter)
            self.density = to_np(rhoElem)
            u, Jelem = self.FE.solve(self.density)

            if (e == 0):
                self.obj0 = ((1.0e-6 + self.density) ** (2 * self.FE.mesh.material['penal']) * Jelem).sum()

            Jelem = (1.0e-6 + self.density) ** (2 * self.FE.mesh.material['penal']) * Jelem
            Jelem = torch.tensor(Jelem).view(-1).float().to(self.device)

            objective = torch.sum(
                torch.div(Jelem, (1.0e-6 + rhoElem) ** self.FE.mesh.material['penal'])) / self.obj0

            volConstraint = ((torch.sum(self.FE.mesh.elemArea * rhoElem) - self.FE.mesh.nelx*self.FE.mesh.nely) / \
                              ((self.FE.mesh.netArea - self.FE.mesh.nelx*self.FE.mesh.nely) * self.desiredVolumeFraction)) - 1.0
            currentVolumeFraction = (torch.sum(self.FE.mesh.elemArea * rhoElem) - self.FE.mesh.nelx*self.FE.mesh.nely) / \
                              (self.FE.mesh.netArea - self.FE.mesh.nelx*self.FE.mesh.nely)
            currentVolumeFraction = currentVolumeFraction.cpu().detach().numpy()

            self.objective = objective
            loss = self.objective + alpha * torch.pow(volConstraint, 2)
            alpha = min(alphaMax, alpha + alphaIncrement)
            loss.backward(retain_graph=True)
            self.optimizer.step()

            compliance.append(self.objective.item() * self.obj0)
            vol.append(currentVolumeFraction)

            if ((e+1) % 30 == 0) and self.densityProjection['sharpness'] < 512:
                if self.densityProjection['sharpness'] < 8:
                    self.densityProjection['sharpness'] = self.densityProjection['sharpness']*2
                else:
                    self.densityProjection['sharpness'] = self.densityProjection['sharpness'] + 4

            self.optimized_latent_code = x
            titleStr = "{:d} \t {:.4F} \t {:.4F}". \
                format(e, self.objective.item() * self.obj0, currentVolumeFraction)
            print(titleStr)

        np.savetxt('convergence_vol.txt', vol)
        np.savetxt('convergence_c.txt', compliance)
        ss = x.cpu().detach().numpy()
        np.savetxt('optimized_latent_code.txt', ss)
        np.savetxt('density.txt', self.density)

if __name__ == "__main__":
    matProp = {'E': 1.0, 'nu': 0.3, 'penal': 3.}  # Structural, starting penal vals
    nelx, nely, nelz = 40, 40, 41
    res = nelx
    elemSize = np.array([1.0, 1.0, 1.0])
    mesh = {'type': 'grid', 'nelx': nelx, 'nely': nely, 'nelz': nelz, 'elemSize': elemSize}
    ndof = 3 * (nelx + 1) * (nely + 1) * (nelz + 1)
    force = np.zeros((ndof, 1))
    fixed = np.array(range(3 * (nely + 1) * (nelx + 1)))
    force_node = np.reshape(range((nely + 1) * (nelx + 1)), (nely + 1, nelx + 1), order='F')
    a = force_node[0::2, 0::2]
    force_node = a.flatten() + (nely + 1) * (nelx + 1) * nelz
    force[3 * force_node + 2, 0] = -1.0 / len(force_node)
    desiredVolumeFraction = 0.30

    bc = {'force': force, 'fixed': fixed}
    densityProjection = {'isOn': True, 'sharpness': 1}
    maxEpochs = 1000

    overrideGPU = True
    start = time.perf_counter()
    topOpt = TopologyOptimizer(mesh, matProp, bc, desiredVolumeFraction, densityProjection, res, overrideGPU)
    topOpt.structural_opt(maxEpochs)

    print("Time taken (secs): {:.2F}".format(time.perf_counter() - start))