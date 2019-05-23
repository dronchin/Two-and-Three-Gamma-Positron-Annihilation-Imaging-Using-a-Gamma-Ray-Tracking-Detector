# -*- coding: utf-8 -*-
"""
Created on Thu May 23 11:08:21 2019

@author: dronchin
This file shows how to open up previously solved three gamma decays.
"""

from ThreeGammaSolver import ThreeGammaImaging


file1 = 'output_good.txt'
file2 = 'output_all_events.txt'

file4 = '3phot_001_output.txt'

file6 = 'ThreePhotCo_001_output.txt'
file7 = 'ThreePhotCo_002_output.txt'
file8 = 'wrong1022Na_001_output.txt'
file9 = 'wrong1022Na_002_output.txt'

file13 = 'wrong1022Na_001_output_with_orig_points.txt'
file14 = 'wrong1022Na_002_output_with_orig_points.txt'
file15 = '3Phot_001_output_with_orig_points.txt'
file16 = 'sim_uniformsource_with_orig_points.txt'
file17 = 'sim_source2_magicpowder_with_orig_points.txt'
file18 = 'testing_file_out.txt'

Reconstructed = ThreeGammaImaging()
Reconstructed.import_previous_run('C:/data/3gamma data sets/{}'.format(file16))
#Reconstructed.solve()

Reconstructed.heatmapplot(20,1, "testfile", False)
