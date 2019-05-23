# -*- coding: utf-8 -*-
"""
Created on Thu May 23 10:48:40 2019

@author: dronchin
Shows how to use ThreeGammaSolver to import a raw data set, solve the decays, then plot a simple heatmap.
"""

from ThreeGammaSolver import ThreeGammaImaging

file1 = '3Phot_001.txt'
file2 = '3Phot_002.txt'
file3 = '3Phot_003.txt'
file4 = '3Phot_004.txt'
file5 = 'ThreePhotFile_00.txt'
file6 = 'ThreePhotFile_00_theta.txt'
file7 = 'sim_pointsource.txt'
file8 = 'sim_source2_magicpowder.txt'



#build the ThreeGammaImaging class
Reconstructed = ThreeGammaImaging()
#import and then solve the three-gamma data
Reconstructed.import_raw_data('C:/data/3gamma data sets/{}'.format(file1))
Reconstructed.solve()
#write solutions to a file
Reconstructed.write_to_file('{}_output.txt'.format(file1[:-4]))

#plot a simple heatmap of the data
Reconstructed.heatmapplot(20,1, "testfile", False)
