#!/usr/bin/env python
# coding: utf-8


from pymatgen.io.vasp import Xdatcar
import numpy as np

from polyhedral_analysis.trajectory import Trajectory
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe, create_matching_site_generator
from polyhedral_analysis.rotation_analyser import RotationAnalyser

import matplotlib.pyplot as plt
from tqdm import tqdm


ps4_recipe = PolyhedraRecipe(central_atoms='P', vertex_atoms='S', method='nearest neighbours',
                            n_neighbours=4)
sns4_recipe = PolyhedraRecipe(central_atoms='Sn', vertex_atoms='S', method='nearest neighbours',
                             n_neighbours=4)

trajectory = Trajectory.from_xdatcars( filenames=['1/XDATCAR',
                                                  '2/XDATCAR',
                                                  '3/XDATCAR',
                                                  '4/XDATCAR',
                                                  '5/XDATCAR'],
                                       recipes=[ps4_recipe, sns4_recipe],
                                       progress=True,
                                       ncores=4 )


reference_points = np.array([[1.0, -1.0, 1.0],
                             [-1.0, -1.0, -1.0],
                             [1.0, 1.0, -1.0],
                             [-1.0, 1.0, 1.0]])
reference_points.shape

all_points = np.array([reference_points, reference_points*(-1)])
all_points.shape

ra = RotationAnalyser(reference_points=all_points)


import math

n = len(trajectory.configurations[0].polyhedra)
# n = 2
all_orientations = []
all_angles = []
all_reference = []
for i in tqdm(range(n)):
    angles = []
    orientations = []
    reference = []
    for c in trajectory.configurations:
        po = ra.polyhedron_orientation( c.polyhedra[i] )
        angles.append(po['rotational_distance'])
        orientations.append(po['orientation_index'])
        reference.append(po['reference_geometry_index'])
    all_orientations.append(orientations)
    all_angles.append(angles)
    all_reference.append(reference)
all_orientations = np.array(all_orientations)
all_angles = np.array(all_angles)*180.0/math.pi



f, ax = plt.subplots(n, 1, sharey=True, figsize=(15,n), gridspec_kw={'hspace': 0})
for axis, angles in zip(ax, all_angles):
    axis.plot(angles)
plt.savefig('angles.png', dpi=300)



f, ax = plt.subplots(n, 1, sharey=True, figsize=(15,n), gridspec_kw={'hspace': 0})
for axis, o in zip(ax, all_orientations):
#     axis.plot(np.not_equal(o[:-1],o[1:]), 'o')
    axis.plot(o)
plt.savefig('orientations.png', dpi=300)



f, ax = plt.subplots(n, 1, sharey=True, figsize=(15,n), gridspec_kw={'hspace': 0})
for axis, o in zip(ax, all_reference):
#     axis.plot(np.not_equal(o[:-1],o[1:]), 'o')
    axis.plot(o)
plt.savefig('reference.png', dpi=300)
