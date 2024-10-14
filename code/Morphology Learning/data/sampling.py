gpu_id = '5'
source_dir = ''
target_dir = ''
################################################

import os
os.environ['CUDA_VISIBLE_DEVICES'] = gpu_id
from datetime import datetime
import numpy as np
import torch
import trimesh
from utils import *
from tqdm import tqdm
from pysdf import SDF

num_samples_and_method = [(500000, 'uniformly'), (500000, 'near')]

num_file = 1
kk = 0
for i in tqdm(range(num_file)):
    t = i
    filename = ''
    in_path = source_dir + filename
    M = trimesh.load_mesh(in_path)
    f = SDF(M.vertices, M.faces)
    mesh = obj2nvc(torch.tensor(M.vertices), torch.tensor(M.faces))
    mesh_normals = face_normals(mesh)
    distrib = area_weighted_distribution(mesh, mesh_normals)

    xyz = sample_points(mesh, num_samples_and_method, mesh_normals, distrib)
    xyz = xyz.cpu().numpy()
    sdf = f(xyz)
    sdf = np.reshape(sdf, (-1, 1))
    xyz_sd = np.concatenate((xyz, sdf), axis=1)
    rand_idx = np.random.permutation(xyz_sd.shape[0])
    xyz_sd = xyz_sd[rand_idx]

    pos = []
    neg = []
    for j in range(xyz_sd.shape[0]):
        if xyz_sd[j, -1] >= 0:
            pos.append(xyz_sd[j, :])
        else:
            neg.append(xyz_sd[j, :])
    outfilename = target_dir + filename + '.npz'
    np.savez(outfilename, pos=np.array(pos), neg=np.array(neg))

