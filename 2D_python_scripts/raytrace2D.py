# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 13:15:33 2022

@author: Rein
"""

import numpy as np
# // Routine:          RayTracing
# // Function:        performs ray tracing to determine attenuation in FOV (not in collimator)
# // Input:                               double& ray tracing angle, double attenuationmap[], int&voxel indices
# // Output:                           Attenuation correction value for input pixel



__all__ = ['RayTracing']

def RayTracing(ThetaRay, AttenuationMap, i, j):
    # // Declare variables
    # double AttenuationCorrectionValue;			// voxel attenuation value ---> output value
    # int delta_ix;								// directional constants used in scan path
    # int delta_iy;								// directional constants used in scan path
    delta_ix = 0
    delta_iy = 0
    
    # // Initial values
    # double s_curr = 0.0;						// current value of s (current position on path)
    # double s_next = 0.0;						// next value of s (next position on path)
    # double LineIntegral = 0.0;					// Line Integral value;
    # double l = 0.0;								// length
    
    s_curr = 0.0						#// current value of s (current position on path)
    s_next = 0.0						#// next value of s (next position on path)
    LineIntegral = 0.0					#// Line Integral value;
    l = 0.0
    
    # // compute unit directions of path
    e_x = np.cos(ThetaRay)
    e_y = np.sin(ThetaRay)
       
    #// prevent dividing by zero
    if e_x == 0:
    	e_x = 1e-15
    if e_y == 0:
    	e_y = 1e-15
    
    #// directional constants used in scan path
    delta_ix = 1 if e_x >= 0 else -1	# determine sign
    d_x = 0 if delta_ix < 0 else 1
    delta_iy = 1 if e_y >= 0 else -1	# determine sign
    d_y = 0 if delta_iy < 0 else 1
    
    #// start voxel
    ix = i
    iy = j
       
    #// distance from start voxel
    Dx = d_x-i
    Dy = d_y-j
    
    #// prevent divisions inside inner loops, so pre-calculate
    inv_e_x = 1/e_x;
    inv_e_y = 1/e_y;
       
    #// compute line integral;
    while (ix > 0) and (iy > 0) and (ix < AttenuationMap.shape[0]-1) and (iy < AttenuationMap.shape[1]-1):
        
        #// possible intersections of path with voxel boundary
        #// (x, y boundary)
        #// direction of the path is taken into account
        s_x = (ix+Dx)*inv_e_x		#//s_x is total path to x-intersection
        s_y = (iy+Dy)*inv_e_y		#//s_y is total path to y-intersection
       
        #// only the closest intersection is really encountered
        #// find this intersection and update voxel index for this direction
        #// in some cases more boundaries are possible (45 degree paths)

        if s_x <= s_y:				#// intersection at x-boundary
        	s_next = s_x
        	ix += delta_ix			#// x-index next voxel
        	
        if s_y <= s_x:				#// intersection at y-boundary
        	s_next = s_y
        	iy += delta_iy			#// y-index next voxel        

        #// length through the voxel
        l = (s_next-s_curr);

        #// calculate line integral
        LineIntegral += AttenuationMap[iy,ix]*l;

        #// update voxelcount and current position
        s_curr = s_next;

    #// end of while loop for ray tracing
   
    #// Calculate attenuation correction value
    return np.exp(-LineIntegral)
 
