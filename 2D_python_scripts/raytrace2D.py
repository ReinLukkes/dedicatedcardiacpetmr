# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 13:15:33 2022

@author: Rein
"""

import numpy as np
import matplotlib.pyplot as plt
import math
# import time
# // Routine:          RayTracing
# // Function:        performs ray tracing to determine attenuation in FOV (not in collimator)
# // Input:                               double& ray tracing angle, double attenuationmap[], int&voxel indices
# // Output:                           Attenuation correction value for input pixel



__all__ = ['RayTracing']

def DetectorLine(radius, detector_num, detector_width, angle, image_size):
    
    detector_array = []
    rad_step = angle - math.pi/2
    center = image_size//2
    
    x0 = center + (np.cos(angle) * radius) - (np.sin(angle) * (detector_width/2))
    y0 = center + (np.sin(angle) * radius) + (np.cos(angle) * (detector_width/2))
    
    detector_array.append( (int(np.round(x0)), int(np.round(y0))) )

    detector_spacing = detector_width / detector_num
    stepx = np.cos(rad_step) * detector_spacing
    stepy = np.sin(rad_step) * detector_spacing


    # returns an array of rounded integer tuples
    for i in range(1, detector_num + 1):
        detector_array.append( ( 
            int(np.round(  
                x0 + (stepx * i))), 
            int(np.round(
                y0 + (stepy * i))) 
                ) )


    return detector_array

def DetectorLineTester(angles, AttenuationMap, detector_num, detector_width, radius, center):
    
    trace_result = np.zeros((detector_num, len(angles) ))
    
    
    print(trace_result.shape[0])
    
    for i, angle in enumerate(np.deg2rad(angles)):
        detector_array = DetectorLine(radius, detector_num, detector_width, -angle, center)
        print(len(detector_array))
        for j in range(detector_num): 
            (x,y) = detector_array[j]

            AttenuationMap[x, y] += 100
            
    
    plt.figure()
    plt.imshow(AttenuationMap, cmap='gray')
    plt.show()


def ForwardProjection2(angles, AttenuationMap, detector_num, detector_width, radius):
    
    trace_result = np.zeros((detector_num, len(angles) ))
    
    
    print(trace_result.shape[0])
    
    for i, angle in enumerate(np.deg2rad(angles)):
        detector_array = DetectorLine(radius, detector_num, detector_width, angle + math.pi, AttenuationMap.shape[0])
        
        # print(angle)
        # print(detector_array)
        
        for j in range(detector_num): 
            (x,y) = detector_array[j]

            trace_result[j, i] = RayTracing(
                angle, 
                AttenuationMap, 
                x, y)
            
            
            # plt.figure()
            # plt.imshow(trace_result, cmap='gray')
            # plt.show()

    return trace_result

def BackProjection2(sinogram, angles, output_size, detector_num, detector_size, detector_width, radius):
    
    trace_result = np.zeros((output_size, output_size))
    DensityMap = np.ones((output_size, output_size))
    
    for i, angle in enumerate(np.deg2rad(angles)):
        detector_array = DetectorLine(radius, detector_num, detector_width, angle + math.pi, output_size)
        
        for j in range(detector_num):
            (x,y) = detector_array[j]
            
            trace_result, DensityMap = inverseRayTracing(
                angle, 
                trace_result, 
                x, y, 
                sinogram[j, i], 
                DensityMap)
    
    # plt.figure()
    # plt.imshow(DensityMap, cmap='gray')
    # plt.show()
    
    return trace_result / DensityMap

def ForwardProjection(angles, AttenuationMap, detector_num, detector_size, detector_height):
    
    trace_result = np.zeros((detector_num, len(angles) ))
    
    print(trace_result.shape[0])
    
    for j in range(detector_num):
        
        for i, angle in enumerate(np.deg2rad(angles)):
            trace_result[j, i] = RayTracing(
                angle, 
                AttenuationMap, 
                AttenuationMap.shape[1]//2 - detector_num//2*detector_size + j*detector_size, 
                detector_height)
    
    # plt.figure()
    # plt.imshow(trace_result, cmap='gray')
    # plt.show()
    
    return trace_result
    
def BackProjection(sinogram, angles, output_size, detector_num, detector_size, detector_height):
    
    trace_result = np.zeros((output_size, output_size))
    DensityMap = np.ones((output_size, output_size))
    
    for j in range(detector_num):
        
        for i, angle in enumerate(np.deg2rad(angles)):
            trace_result, DensityMap = inverseRayTracing(
                angle, 
                trace_result, 
                trace_result.shape[1]//2 - detector_num//2*detector_size + j*detector_size, 
                detector_height, 
                sinogram[j, i], 
                DensityMap)
    
    # plt.figure()
    # plt.imshow(DensityMap, cmap='gray')
    # plt.show()
    
    return trace_result / DensityMap    


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

    # print("Angle =", ThetaRay)
    #// compute line integral;
    while (ix > 0) and (iy > 0) and (ix < AttenuationMap.shape[1]-1) and (iy < AttenuationMap.shape[0]-1):
        
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
        l = (s_next-s_curr)

        #// calculate line integral
        LineIntegral += AttenuationMap[iy,ix]*l

        #// update voxelcount and current position
        s_curr = s_next

    #// end of while loop for ray tracing
        

    #// Calculate attenuation correction value
    # if LineIntegral != 0 and np.exp(-LineIntegral) != 0:
    #     print("Set:")
    #     print(LineIntegral)
    #     print(np.exp(-LineIntegral))
    #     time.sleep(.1)
    return LineIntegral
    # return np.exp(-LineIntegral) 


def inverseRayTracing(ThetaRay, AttenuationMap, i, j, value, DensityMap):
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
    
    # print("Angle =", ThetaRay)
    #// compute line integral;
    while (ix > 0) and (iy > 0) and (ix < AttenuationMap.shape[1]-1) and (iy < AttenuationMap.shape[0]-1):
        
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
        l = (s_next-s_curr)

        #// update the value of the pixel  
        
        AttenuationMap[iy, ix] += value
        
        # add value to density map
        DensityMap[iy, ix] += .1

        #// update voxelcount and current position
        s_curr = s_next

    #// end of while loop for ray tracing

    return AttenuationMap, DensityMap