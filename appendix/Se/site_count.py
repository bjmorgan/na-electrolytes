#!/usr/bin/env python
# coding: utf-8

from pymatgen.io.vasp import Poscar, Xdatcar
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen import Structure, Lattice
import numpy as np
import operator
from site_analysis import Atom, Analysis, ShortestDistanceSite, get_vertex_indices, AtomsTrajectory, SitesTrajectory
from collections import Counter
import tqdm

x1 = Xdatcar('1/XDATCAR')
x2 = Xdatcar('2/XDATCAR')
x3 = Xdatcar('3/XDATCAR')
x4 = Xdatcar('4/XDATCAR')
x5 = Xdatcar('5/XDATCAR')
structures = x1.structures + x2.structures + x3.structures + x4.structures + x5.structures


all_na_structure = Poscar.from_file('na_sn_all_na_ext.POSCAR.vasp').structure
vertex_species = 'Se'
centre_species = 'Na'

sg = SpaceGroup('I41/acd:2')
from pymatgen import Structure, Lattice
lattice = all_na_structure.lattice
na1 = Structure.from_spacegroup(sg='I41/acd:2', lattice=lattice, species=['Na'], coords=[[0.25, 0.0, 0.125]])
na2 = Structure.from_spacegroup(sg='I41/acd:2', lattice=lattice, species=['Na'], coords=[[0.00, 0.0, 0.125]])
na3 = Structure.from_spacegroup(sg='I41/acd:2', lattice=lattice, species=['Na'], coords=[[0.0, 0.25, 0.0]])
na4 = Structure.from_spacegroup(sg='I41/acd:2', lattice=lattice, species=['Na'], coords=[[0.0, 0.0, 0.0]])
na5 = Structure.from_spacegroup(sg='I41/acd:2', lattice=lattice, species=['Na'], coords=[[0.75, 0.25, 0.0]])
na6 = Structure.from_spacegroup(sg='I41/acd:2', lattice=lattice, species=['Na'], coords=[[0.5, 0.75, 0.625]])
i2  = Structure.from_spacegroup(sg='I41/acd:2', lattice=lattice, species=['Na'], coords=[[0.666, 0.1376, 0.05]])
na_structures = {'Na1': na1,
                 'Na2': na2,
                 'Na3': na3,
                 'Na4': na4,
                 'Na5': na5,
                 'Na6': na6,
                 'i2': i2}

na1_sites = [ ShortestDistanceSite(s.frac_coords, label='Na1') for s in na1 ]
na2_sites = [ ShortestDistanceSite(s.frac_coords, label='Na2') for s in na2 ]
na3_sites = [ ShortestDistanceSite(s.frac_coords, label='Na3') for s in na3 ]
na4_sites = [ ShortestDistanceSite(s.frac_coords, label='Na4') for s in na4 ]
na5_sites = [ ShortestDistanceSite(s.frac_coords, label='Na5') for s in na5 ]
na6_sites = [ ShortestDistanceSite(s.frac_coords, label='Na6') for s in na6 ]
i2_sites  = [ ShortestDistanceSite(s.frac_coords, label='i2') for s in i2 ]
sites = na1_sites + na2_sites + na3_sites + na4_sites + na5_sites + na6_sites + i2_sites


structure = Poscar.from_file('1/POSCAR').structure
# create Polyhedron objects
# create Atom objects
atoms = [Atom(species_string=centre_species) for site in structure if site.species_string is 'Na']
analysis = Analysis(sites, atoms)

analysis.trajectory_from_structures( structures, progress=True)


#Counts instances of no, single, and double occupation for each site.
n_timesteps = len(analysis.timesteps)
c_sites = { l: Counter() for l in analysis.site_labels() }
c = Counter()
p_occ = {}
for site in analysis.sites:
    for ts in site.trajectory:
        c_sites[site.label][len(ts)] += 1
f = open("detail_sites.dat", "w+")
f.write(str(c_sites))
f.close()


#Altered version of summation for probabilities that includes sum of len in order to account for double site occupations. This was causing the inital bug.
n_timesteps = len(analysis.timesteps)
c_sites = Counter(analysis.site_labels())
c = Counter()
p_occ = {}
for site in analysis.sites:
    c[site.label] += sum([ len(ts) for ts in site.trajectory if len(ts)>0 ])

for k, v in c.items():
    p_occ[k] = v / c_sites[k] / n_timesteps
p_occ
f = open("sites.dat", "w+")
f.write(str(p_occ))
f.close()


for k,v in c.items():
    check = sum( [ p_occ[k] * c_sites[k] for k, v in c.items()])
f = open("check_sites.dat", "w+")
f.write(str(check))
f.close()
