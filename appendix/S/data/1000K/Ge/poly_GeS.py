#!/bin/python

from pymatgen.io.vasp import Xdatcar
import numpy as np

from polyhedral_analysis.trajectory import Trajectory
from polyhedral_analysis.polyhedra_recipe import PolyhedraRecipe, create_matching_site_generator
from polyhedral_analysis.rotation_analyser import RotationAnalyser

import matplotlib.pyplot as plt
from tqdm import tqdm

ps4_recipe = PolyhedraRecipe(central_atoms='P', vertex_atoms='S', method='nearest neighbours',
                            n_neighbours=4)
sns4_recipe = PolyhedraRecipe(central_atoms='Ge', vertex_atoms='S', method='nearest neighbours',
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

inv_all_points = np.array([reference_points*(-1), reference_points])
inv_all_points.shape

ra = RotationAnalyser(reference_points=all_points)
inv_ra = RotationAnalyser(reference_points=inv_all_points)

import math

type_0 = [2,3,4,5]
all_orientations = []
all_angles = []
all_reference = []

type_1 = [0,1,6,7]
inv_all_orientations = []
inv_all_angles = []
inv_all_reference = []

for i in tqdm(type_0):
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
all_reference = np.array(all_reference)

for i in tqdm(type_1):
    angles = []
    orientations = []
    reference = []
    for c in trajectory.configurations:
        po = inv_ra.polyhedron_orientation( c.polyhedra[i] )
        angles.append(po['rotational_distance'])
        orientations.append(po['orientation_index'])
        reference.append(po['reference_geometry_index'])
    inv_all_orientations.append(orientations)
    inv_all_angles.append(angles)
    inv_all_reference.append(reference)
inv_all_orientations = np.array(inv_all_orientations)
inv_all_angles = np.array(inv_all_angles)*180.0/math.pi
inv_all_reference = np.array(inv_all_reference)


from skimage.util import view_as_windows

def smooth_data(a, n):
    return np.array([ np.argmax(np.bincount(window)) for window in view_as_windows(a,n) ])

#reference graphs
n = 4
smoothed_data = [smooth_data(all_reference[i], 200) for i in range(4)]
sosmooth_data = [smooth_data(smoothed_data[i], 100) for i in range(4)]
f, ax = plt.subplots(n, 1, sharey=True, figsize=(15,n), gridspec_kw={'hspace': 0})
for axis, o in zip(ax, sosmooth_data):
    axis.plot(o)
ax[0].set_title('Reference State of α-PS4 Tetrahedra')
ax[3].set_xlabel('Time / fs')
f.text(0.1, 0.5, 'Ground State = 0, Excited State = 1', va='center', rotation='vertical')
for i in range(3):
    ax[i].set_xticks([])
plt.savefig('reference_0.png', dpi=300)

n = 4
smoothed_data = [smooth_data(inv_all_reference[i], 200) for i in range(4)]
sosmooth_data = [smooth_data(smoothed_data[i], 100) for i in range(4)]
f, ax = plt.subplots(n, 1, sharey=True, figsize=(15,n), gridspec_kw={'hspace': 0})
for axis, o in zip(ax, sosmooth_data):
    axis.plot(o)
ax[0].set_title('Reference State of β-PS4 Tetrahedra')
ax[3].set_xlabel('Time / fs')
f.text(0.1, 0.5, 'Ground State = 0, Excited State = 1', va='center', rotation='vertical')
for i in range(3):
    ax[i].set_xticks([])
plt.savefig('reference_1.png', dpi=300)

#orientation graphs
n = 4
smoothed_data = [smooth_data(all_orientations[i], 200) for i in range(4)]
sosmooth_data = [smooth_data(smoothed_data[i], 100) for i in range(4)]
f, ax = plt.subplots(n, 1, sharey=True, figsize=(15,n), gridspec_kw={'hspace': 0})
for axis, o in zip(ax, sosmooth_data):
    axis.plot(o)
ax[0].set_title('Orientation of α-PS4 Tetrahedra')
ax[3].set_xlabel('Time / fs')
f.text(0.1, 0.5, 'Closest Orientation Number', va='center', rotation='vertical')
for i in range(3):
    ax[i].set_xticks([])
plt.savefig('orientations_0.png', dpi=300)

n = 4
smoothed_data = [smooth_data(inv_all_orientations[i], 200) for i in range(4)]
sosmooth_data = [smooth_data(smoothed_data[i], 100) for i in range(4)]
f, ax = plt.subplots(n, 1, sharey=True, figsize=(15,n), gridspec_kw={'hspace': 0})
for axis, o in zip(ax, sosmooth_data):
    axis.plot(o)
ax[0].set_title('Orientation of β-PS4 Tetrahedra')
ax[3].set_xlabel('Time / fs')
f.text(0.1, 0.5, 'Closest Orientation Number', va='center', rotation='vertical')
for i in range(3):
    ax[i].set_xticks([])
plt.savefig('orientations_1.png', dpi=300)

#angle graphs
n = 4
f, ax = plt.subplots(n, 1, sharey=True, figsize=(15,n), gridspec_kw={'hspace': 0})
for axis, o in zip(ax, all_angles):
    axis.plot(o)
ax[0].set_title('Vibration of α-PS4 Tetrahedra')
ax[3].set_xlabel('Time / fs')
f.text(0.09, 0.5, 'Angle to Closest Orientation / °', va='center', rotation='vertical')
for i in range(3):
    ax[i].set_xticks([])
plt.savefig('angle_0.png', dpi=300)

n = 4
f, ax = plt.subplots(n, 1, sharey=True, figsize=(15,n), gridspec_kw={'hspace': 0})
for axis, o in zip(ax, inv_all_angles):
    axis.plot(o)
ax[0].set_title('Vibration of β-PS4 Tetrahedra')
ax[3].set_xlabel('Time / fs')
f.text(0.09, 0.5, 'Angle to Closest Orientation / °', va='center', rotation='vertical')
for i in range(3):
    ax[i].set_xticks([])
plt.savefig('angle_1.png', dpi=300)


#reference data
type_0 = [0,1,2,3]
type_0_rotation_count = 0
type_0_poly = 0
for n in type_0:
    smooth = smooth_data(all_reference[n], 200)
    sosmooth = smooth_data(smooth, 100)
    for j, k in enumerate(sosmooth):
        if j == len(sosmooth) - 1:
            type_0_poly += 1
        elif j == 0:
            pass
        elif k == sosmooth[j+1] and k != sosmooth[j-1]:
            type_0_rotation_count += 1
with open('reference_0.dat', 'w+') as f:
    f.write('α-PS4 State Changes' + ' ' + str(type_0_rotation_count) + '\n'
            + 'α-PS4 Tetrahedra' + ' ' + str(type_0_poly) + '\n'
            + 'Mean α-PS4 State Change Frequency' + ' ' + str(type_0_rotation_count/type_0_poly) + ' '
            + 'fs\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}')

type_1 = [0,1,2,3]
type_1_rotation_count = 0
type_1_poly = 0
for n in type_0:
    smooth = smooth_data(inv_all_reference[n], 200)
    sosmooth = smooth_data(smooth, 100)
    for j, k in enumerate(sosmooth):
        if j == len(sosmooth) - 1:
            type_1_poly += 1
        elif j == 0:
            pass
        elif k == sosmooth[j+1] and k != sosmooth[j-1]:
            type_1_rotation_count += 1
with open('reference_1.dat', 'w+') as f:
    f.write('β-PS4 State Changes' + ' ' + str(type_1_rotation_count) + '\n'
            + 'β-PS4 Tetrahedra' + ' ' + str(type_1_poly) + '\n'
            + 'Mean β-PS4 State Change Frequency' + ' ' + str(type_1_rotation_count/type_1_poly) + ' '
            + 'fs\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}')

#orientation data
type_0 = [0,1,2,3]
type_0_rotation_count = 0
type_0_poly = 0
for n in type_0:
    smooth = smooth_data(all_orientations[n], 200)
    sosmooth = smooth_data(smooth, 100)
    for j, k in enumerate(sosmooth):
        if j == len(sosmooth) - 1:
            type_0_poly += 1
        elif j == 0:
            pass
        elif k == sosmooth[j+1] and k != sosmooth[j-1]:
            type_0_rotation_count += 1
with open('orientation_0.dat', 'w+') as f:
    f.write('α-PS4 Orientation Changes' + ' ' + str(type_0_rotation_count) + '\n'
            + 'α-PS4 Tetrahedra' + ' ' + str(type_0_poly) + '\n'
            + 'Mean α-PS4 Rotation Frequency' + ' ' + str(type_0_rotation_count/type_0_poly) + ' '
            + 'fs\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}')

type_1 = [0,1,2,3]
type_1_rotation_count = 0
type_1_poly = 0
for n in type_0:
    smooth = smooth_data(inv_all_orientations[n], 200)
    sosmooth = smooth_data(smooth, 100)
    for j, k in enumerate(sosmooth):
        if j == len(sosmooth) - 1:
            type_1_poly += 1
        elif j == 0:
            pass
        elif k == sosmooth[j+1] and k != sosmooth[j-1]:
            type_1_rotation_count += 1
with open('orientation_1.dat', 'w+') as f:
    f.write('β-PS4 Orientation Changes' + ' ' + str(type_1_rotation_count) + '\n'
            + 'β-PS4 Tetrahedra' + ' ' + str(type_1_poly) + '\n'
            + 'Mean Type β-PS4 Rotation Frequency' + ' ' + str(type_1_rotation_count/type_1_poly) + ' '
            + 'fs\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}')

#ground state vs excited state
type_0 = [0,1,2,3]
time_0_0 = 0
time_0_1 = 0
for n in type_0:
    for j, k in enumerate(all_reference[n]):
        if k == 0:
            time_0_0 += 1
        elif k == 1:
            time_0_1 += 1

type_1 = [0,1,2,3]
time_1_0 = 0
time_1_1 = 0
for n in type_1:
    for j, k in enumerate(inv_all_reference[n]):
        if k == 0:
            time_1_0 += 1
        elif k == 1:
            time_1_1 += 1

total_time = time_0_0 + time_1_0 + time_0_1 + time_1_1
ground_percent = (time_0_0 + time_1_0) / total_time * 100
excited_percent = (time_0_1 + time_1_1) / total_time * 100

with open('ground_excited.dat', 'w+') as f:
    f.write('α-PS4 Time in ground state' + ' ' + str(time_0_0 * 2) + ' ' + 'fs' + '\n'
            + 'α-PS4 Time in excited state' + ' ' + str(time_0_1 * 2) + ' ' + 'fs' + '\n'
            + 'β-PS4 Time in ground state' + ' ' + str(time_1_0 * 2) + ' ' + 'fs' + '\n'
            + 'β-PS4 Time in excited state' + ' ' + str(time_1_1 * 2) + ' ' + 'fs' + '\n'
            + 'PS4 Percentage time in ground state' + ' ' + str(ground_percent) + ' ' + '%' + '\n'
            + 'PS4 Percentage time in excited state' + ' ' + str(excited_percent) + ' ' + '%' + '\n')
