#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
from pymatgen import Structure
from pymatgen.io.vasp import Poscar, Xdatcar
from pymatgen.analysis.diffusion_analyzer import DiffusionAnalyzer, get_arrhenius_plot, get_extrapolated_conductivity

temperatures = [500, 600, 700, 800, 900, 1000]
tetrels = ["Ge", "Sn"]

for i in temperatures:
    for j in tetrels:
        x1 = Xdatcar(f'{i}K/{j}/1/XDATCAR')
        x2 = Xdatcar(f'{i}K/{j}/2/XDATCAR')
        x3 = Xdatcar(f'{i}K/{j}/3/XDATCAR')
        x4 = Xdatcar(f'{i}K/{j}/4/XDATCAR')
        x5 = Xdatcar(f'{i}K/{j}/5/XDATCAR')
        structures = x1.structures + x2.structures + x3.structures + x4.structures + x5.structures
        


s500K = ["500/Sn/1/vasprun.xml", "500/Sn/2/vasprun.xml", "500/Sn/3/vasprun.xml"]
d500K = DiffusionAnalyzer.from_files(s500K, specie="Na", smoothed=max)
s600K = ["600/Sn/1/vasprun.xml", "600/Sn/2/vasprun.xml", "600/Sn/3/vasprun.xml"]
d600K = DiffusionAnalyzer.from_files(s600K, specie="Na", smoothed=max)
s700K = ["700/Sn/1/vasprun.xml", "700/Sn/2/vasprun.xml", "700/Sn/3/vasprun.xml"]
d700K = DiffusionAnalyzer.from_files(s700K, specie="Na", smoothed=max)
s800K = ["800/Sn/1/vasprun.xml", "800/Sn/2/vasprun.xml", "800/Sn/3/vasprun.xml"]
d800K = DiffusionAnalyzer.from_files(s800K, specie="Na", smoothed=max)
s900K = ["900/Sn/1/vasprun.xml", "900/Sn/2/vasprun.xml", "900/Sn/3/vasprun.xml"]
d900K = DiffusionAnalyzer.from_files(s900K, specie="Na", smoothed=max)
s1000K = ["1000/Sn/1/vasprun.xml", "1000/Sn/2/vasprun.xml", "1000/Sn/3/vasprun.xml"]
d1000K = DiffusionAnalyzer.from_files(s1000K, specie="Na", smoothed=max)

temperatures = [500, 600, 700, 800, 900, 1000]
analyzers = [d500K, d600K, d700K, d800K, d900K, d1000K]
diffusivities = [d.diffusivity for d in analyzers]
plt = get_arrhenius_plot(temperatures, diffusivities)
plt.savefig('sn_na_arrhenius_plot.png', dpi=300)

rts = get_extrapolated_conductivity(temperatures, diffusivities, 
                                    new_temp=300, structure=analyzers[0].structure, 
                                    species="Na")
print("The Na ionic conductivity for Na11Sn2PS12 at 300 K is %.4f mS/cm" % rts)
