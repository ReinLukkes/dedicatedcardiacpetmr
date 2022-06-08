## Generate lookup table PET-detectors in UMCU PET-MRI system
# Note: all lengths are in cm 

#       Y                                        _________  
#       |                                       / _ \     \ 
#       |                                      | / \ |     |
#       |_____ Z                               | | | |     |
#        \                                     | | | |     |
#         \                                    | \_/ |     |
#          X                                    \___/_____/

## Information PET-system
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import customLORindices

n_sectors = 12             # number of sectors
R = 38.5 + 1.6/2           # radius PET-system
n_mod_ax = 3               # number of modules axial
n_mod_transax = 3          # number of modules transaxial
gap_mod_ax = 0.207         # module gap axial
gap_mod_transax = 0.807    # module gap transaxial
n_det_per_mod_ax = 12      # number of detectors per module axial
n_det_per_mod_transax = 12 # number of detectors per module transaxial
det_size = 0.393           # size of detector axial and transaxial
gap_det = 0.007            # detector gap axial and transaxial

## 

n_rings = n_mod_ax * n_det_per_mod_ax                          # number of rings
n_det_per_sect_transax = n_mod_transax * n_det_per_mod_transax # number of detectors per sector in one ring
n_det_per_ring = n_det_per_sect_transax * n_sectors            # number of detectors per ring
d_angle = 360.0 / n_sectors                                    # step size of angle to middle of a sector

PET_PHILIPS_UMCU_geom = np.zeros((n_rings,n_det_per_ring,4))                 # lookup table of the geometry

z = .5 * det_size  # z of first ring of detectors

# loop over all rings
# for r = 1:n_rings:
for r in range(n_rings):
    
    # loop over all sectors
    # for i = 1:n_sectors:
    for i in range(n_sectors):
        
        angle = np.deg2rad(90 + 8 + (i) * d_angle)    # angle sector, verified in GATE
        angle_sector = angle + (np.pi/2)      # slope of sector
        x_mid = R * np.cos(angle)        # x middle of sector
        y_mid = R * np.sin(angle)        # y middle of sector
        
        x = x_mid + .5 * (det_size + gap_det) * np.cos(angle_sector) # x of first detector after middle
        y = y_mid + .5 * (det_size + gap_det) * np.sin(angle_sector) # y of first detector after middle
               
        det = (i + .5) * n_det_per_sect_transax + 1        # detector index first detector after middle
        
        # loop over all detectors after middle of current sector
        while det < ((i+1) * n_det_per_sect_transax + 1):
            
            # add x, y and z to lookup table for current detector
            PET_PHILIPS_UMCU_geom[r,det,1] = x
            PET_PHILIPS_UMCU_geom[r,det,2] = y
            PET_PHILIPS_UMCU_geom[r,det,3] = z
            PET_PHILIPS_UMCU_geom[r,det,4] = (r)* n_det_per_ring + det
            
            # check if the next detector is in a different module
            if round(det / n_det_per_mod_transax) == (det / n_det_per_mod_transax):
                x = x + (det_size + gap_mod_transax) * np.cos(angle_sector)  # add module gap
                y = y + (det_size + gap_mod_transax) * np.sin(angle_sector)  # add module gap
            else:
                x = x + (det_size + gap_det) * np.cos(angle_sector)  # only add detector gap
                y = y + (det_size + gap_det) * np.sin(angle_sector)  # only add detector gap
            
            det = det + 1  # index next detector
        
        x = x_mid - .5 * (det_size + gap_det) * np.cos(angle_sector) # x of first detector before middle
        y = y_mid - .5 * (det_size + gap_det) * np.sin(angle_sector) # y of first detector before middle
        
        det = (i + .5) * n_det_per_sect_transax             # detector index first detector before middle
        
        # loop over all detectors before middle of current sector
        while det > (i) * n_det_per_sect_transax:
            
            # add x, y and z to lookup table for current detector
            PET_PHILIPS_UMCU_geom[r,det,1] = x
            PET_PHILIPS_UMCU_geom[r,det,2] = y
            PET_PHILIPS_UMCU_geom[r,det,3] = z
            PET_PHILIPS_UMCU_geom[r,det,4] = (r)*n_det_per_ring + det
            
            # check if the next detector is in a different module
            if round(((det - 1) / n_det_per_mod_transax)) == ((det - 1) / n_det_per_mod_transax):
                x = x - (det_size + gap_mod_transax) * np.cos(angle_sector)  # add module gap
                y = y - (det_size + gap_mod_transax) * np.sin(angle_sector)  # add module gap
            else:
                x = x - (det_size + gap_det) * np.cos(angle_sector)  # only add detector gap
                y = y - (det_size + gap_det) * np.sin(angle_sector)  # only add detector gap
            
            det = det - 1  # index next detector
        
    
    # check if the next ring of detectors is in a different module
    if ((r+1) // n_det_per_mod_ax) == 1 or ((r+1) // n_det_per_mod_ax) == 2:
        z = z + det_size + gap_mod_ax  # add module gap
    else:
        z = z + det_size + gap_det # only add detector gap
    

# PET_PHILIPS_UMCU_geom = circshift(PET_PHILIPS_UMCU_geom, -1*n_det_per_mod_ax*n_mod_ax, 2)
# PET_PHILIPS_UMCU_geom = flip(PET_PHILIPS_UMCU_geom,2)
# PET_PHILIPS_UMCU_geom = flip(PET_PHILIPS_UMCU_geom,1)
# [LOR_indices, ~] = getLORindices(PET_PHILIPS_UMCU_geom)


LOR_indices = customLORindices(PET_PHILIPS_UMCU_geom, n_sectors, n_mod_ax, n_mod_transax, n_det_per_mod_transax, n_det_per_mod_ax)

datafile = pd.DataFrame(LOR_indices)
datafile.to_csv('LOR_indices.csv', header = False, index = False)

#save('C:\Users\rjosesan\Documents\MATLAB\SSS_Algorithm\ValidatedGeom\geom.mat', 'PET_PHILIPS_UMCU_geom')
#save('C:\Users\rjosesan\Documents\MATLAB\SSS_Algorithm\ValidatedGeom\lor.mat', 'LOR_indices')