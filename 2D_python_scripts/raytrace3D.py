# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:52:41 2022

@author: rlukkes2
"""

import numpy as np
import matplotlib.pyplot as plt
import math
# import time
# // Routine:          RayTracing
# // Function:        performs ray tracing to determine attenuation in FOV (not in collimator)
# // Input:                               double& ray tracing angle, double attenuationmap[], int&voxel indices
# // Output:                           Attenuation correction value for input pixel

# Value added to the density map on each pass through a pixel. Used to normalise results
densityMapStrength = 0.1

class Detector:
    
    def __init__(self, radius, number, detectorType, angles, attenuationMap, width = 0, cellSize = 1):
        self.radius     = radius    # Distance of scanner from centre
        self.number     = number    # Number of cells in scanner
        self.cellSize   = cellSize  # Size of cells in dedicated scanner
        self.detectorType      = detectorType     # Detector type: 0 = dedicated, 1 = full
        self.angles     = angles
        self.attenuationMap = attenuationMap
        if width == 0:              # width of the partial detector array used for full
            self.width = number
        else:
            self.width = width
            
        if detectorType != 0 and detectorType != 1:
            print("error: detectorType invalid")
            
        if number <= 0:
            print("error: detector number too small")
            
        if detectorType == 1:
            self.detectorArray = self.CircleDetector()
            
            
    def sfp(self):  
        if self.detectorType == 0:
            return self.ForwardProjection()
        if self.detectorType == 1:
            return self.ForwardProjection2()
            
        
    def sbp(self, sinogram):
        if self.detectorType == 0:
            return self.BackProjection(sinogram)
        if self.detectorType == 1:
            return self.BackProjection2(sinogram)

    def setAttenuationMap(self, am):
        self.attenuationMap = am
        
    __all__ = ['sfp', 'sbp', 'setAttenuationMap']
    
    def CircleDetector(self):
        
        detector_array = []
        angles = np.arange(0, 2*np.pi, 2*np.pi/self.number)
        center = self.attenuationMap.shape[0]//2
        

        for theta in angles:
            x0 = center + ( (np.cos(theta) * self.radius))
            y0 = center + ( (np.sin(theta) * self.radius))
            detector_array.append(
                (int(np.round(x0)), int(np.round(y0)), 0)
                )
    
        return detector_array
    
    
    def ForwardProjection2(self):
        
        sinogram = np.zeros((self.number, self.number))
        
        for i in range(self.number):
            if i%25 == 0:
                print(i, '/', self.number)
            for j in range(self.number):
                if i == j:
                    continue
                sinogram[j, i] = self.RayTracing(j, i)
                
        # plt.figure()
        # plt.imshow(sinogram, cmap='gray')
        # plt.show()        
        
        return sinogram
    
    def BackProjection2(self, sinogram):
        dim = self.attenuationMap.shape[0]
        z = self.attenuationMap.shape[2]
        trace_result = np.zeros((dim, dim, z))
        DensityMap = np.ones((dim, dim, z))


        for i in range(self.number):
            if i%25 == 0:
                print(i, '/', self.number)
            for j in range(self.number):
                if i == j:
                    continue
                
                trace_result, DensityMap = self.inverseRayTracing(
                    trace_result, 
                    DensityMap, 
                    j, i, 
                    sinogram[j, i]
                    )
                
        # plt.figure()
        # plt.imshow(DensityMap, cmap='gray')
        # plt.show()
        
        return trace_result / DensityMap
    
    def ForwardProjection(self, AttenuationMap):
        
        trace_result = np.zeros((self.number, len(self.angles) ))
        
        #print(trace_result.shape[0])
        
        for j in range(self.number):
            
            for i, angle in enumerate(np.deg2rad(self.angles)):
                trace_result[j, i] = self.RayTracing(
                    angle, 
                    AttenuationMap, 
                    AttenuationMap.shape[1]//2 - self.number//2*self.cellSize + j*self.cellSize, 
                    AttenuationMap.shape[1]//2 - self.radius)
        
        # plt.figure()
        # plt.imshow(trace_result, cmap='gray')
        # plt.show()
        
        return trace_result
        
    def BackProjection(self, sinogram, output_size):
        
        trace_result = np.zeros((output_size, output_size))
        DensityMap = np.ones((output_size, output_size))
        
        for j in range(self.number):
            
            for i, angle in enumerate(np.deg2rad(self.angles)):
                trace_result, DensityMap = self.inverseRayTracing(
                    angle, 
                    trace_result, 
                    output_size//2 - self.number//2*self.cellSize + j*self.cellSize, 
                    output_size//2 - self.radius, 
                    sinogram[j, i], 
                    DensityMap)
        
        # plt.figure()
        # plt.imshow(DensityMap, cmap='gray')
        # plt.show()
        
        return trace_result / DensityMap    
    
    
    def RayTracing(self, j, i):
        # // Declare variables
        # double AttenuationCorrectionValue;			// voxel attenuation value ---> output value
        # int delta_ix;								// directional constants used in scan path
        # int delta_iy;								// directional constants used in scan path
        delta_ix = 0
        delta_iy = 0
        delta_iz = 0
        
        startVoxel = self.detectorArray[j]
        endVoxel = self.detectorArray[i]
        rayVector = tuple(map(lambda k, l: k-l, endVoxel, startVoxel))      # get vector between the two voxels
        magnitude = np.sqrt(np.sum(list(map(lambda k: k**2, rayVector))))   # calculate the magnitude
        rayVector = tuple(map(lambda k: k/magnitude, rayVector))            # calculate the unit vector
        rayVector = tuple(map(lambda k: 1e-15 if k == 0 else k, rayVector))            # prevent dividing by zero

        mapDim = self.attenuationMap.shape

        # // Initial values
        s_curr = 0.0						#// current value of s (current position on path)
        s_next = 0.0						#// next value of s (next position on path)
        LineIntegral = 0.0					#// Line Integral value;
        l = 0.0        
        
        # TEMPORARY
        R_FOV = 1
        
        # % assuming isotropic voxels for calculation of path length (in cm)
        # voxel_dim = R_FOV / (.5 * mapDim[1]);
        
        voxel_dim = 1
        
        # set voxel to start
        ix = startVoxel[0]
        iy = startVoxel[1]
        iz = startVoxel[2]
        
        #// directional constants used in scan path
        delta_ix = 1 if rayVector[0] >= 0 else -1	# determine sign
        d_x = 0 if delta_ix < 0 else 1
        delta_iy = 1 if rayVector[1] >= 0 else -1	# determine sign
        d_y = 0 if delta_iy < 0 else 1
        delta_iz = 1 if rayVector[2] >= 0 else -1	# determine sign
        d_z = 0 if delta_iz < 0 else 1
        
        #// distance from start voxel
        Dx = d_x-ix
        Dy = d_y-iy
        Dz = d_z-iz
        
        #// prevent divisions inside inner loops, so pre-calculate
        inv_e_x = 1/rayVector[0]
        inv_e_y = 1/rayVector[1]
        inv_e_z = 1/rayVector[2]
    
        #// compute line integral;
        while (ix >= 0) and (iy >= 0) and (iz >= 0) and (ix < mapDim[1]-1) and (iy < mapDim[0]-1) and (iz < mapDim[2]-1):
            
            #// possible intersections of path with voxel boundary
            #// (x, y boundary)
            #// direction of the path is taken into account
            s_x = (ix+Dx)*inv_e_x		#//s_x is total path to x-intersection
            s_y = (iy+Dy)*inv_e_y		#//s_y is total path to y-intersection
            s_z = (iz+Dz)*inv_e_z		#//s_z is total path to z-intersection
           
            # % calculate distance to each voxel boundary
            # a_x = abs((ix-x)*inv_e_x);
            # a_y = abs((iy-y)*inv_e_y);
            # a_z = abs((iz-z)*inv_e_z);
            
            
            #// only the closest intersection is really encountered
            #// find this intersection and update voxel index for this direction
            #// in some cases more boundaries are possible (45 degree paths)
            #TDOD: check influence of x voxel preference
            
            
            
            if s_x <= s_y and s_x <= s_z:				#// intersection at x-boundary
               	s_next = s_x
               	ix += delta_ix			#// x-index next voxel
                if s_x == s_y: iy += delta_iy 
                if s_x == s_z: iz += delta_iz

            elif s_y <= s_x and s_y <= s_z:				#// intersection at y-boundary
                s_next = s_y
                iy += delta_iy			#// y-index next voxel
                if s_y == s_z: iz += delta_iz

            else: #s_z <= s_x and s_z <= s_y
                s_next = s_z
                iz += delta_iz
    
    
    
            #TODO: proper length through voxel calculation
            
            
            # % calculate the length through current voxel
            # l = (((1+(e_z^2))*(a^2))^.5)*voxel_dim;
            l = np.sqrt( (1 + (rayVector[2]**2)) * (s_next**2) ) * voxel_dim
    
            # % get the attenuation coefficient
            # if AttenuationMap(round(x),round(y),round(z)) > 0
            #     tissue = AttenuationMap(round(x),round(y),round(z));
            #     atn_coef = attenuation_table(energy_index,tissue);
            # else
            #     atn_coef = 0;
            # end
            
            # LineIntegralAtt = LineIntegralAtt + l*atn_coef;
            # LineIntegralWeight = LineIntegralWeight + l*ActivityMap(round(x),round(y),round(z));
    
            #// calculate line integral
            LineIntegral += self.attenuationMap[iy,ix,iz]*l 
    
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
    
    
    def inverseRayTracing(self, AttenuationMap, DensityMap, j, i, sinogramVal):
        
        # // Declare variables
        # double AttenuationCorrectionValue;			// voxel attenuation value ---> output value
        # int delta_ix;								// directional constants used in scan path
        # int delta_iy;								// directional constants used in scan path
        delta_ix = 0
        delta_iy = 0
        delta_iz = 0
        
        startVoxel = self.detectorArray[j]
        endVoxel = self.detectorArray[i]
        rayVector = tuple(map(lambda k, l: k-l, endVoxel, startVoxel))      # get vector between the two voxels
        magnitude = np.sqrt(np.sum(list(map(lambda k: k**2, rayVector))))   # calculate the magnitude
        rayVector = tuple(map(lambda k: k/magnitude, rayVector))            # calculate the unit vector
        rayVector = tuple(map(lambda k: 1e-15 if k == 0 else k, rayVector))            # prevent dividing by zero

        mapDim = self.attenuationMap.shape

        # // Initial values
        s_curr = 0.0						#// current value of s (current position on path)
        s_next = 0.0						#// next value of s (next position on path)
        LineIntegral = 0.0					#// Line Integral value;
        l = 0.0        
        
        # TEMPORARY
        R_FOV = 1
        
        # % assuming isotropic voxels for calculation of path length (in cm)
        # voxel_dim = R_FOV / (.5 * mapDim[1]);
        
        voxel_dim = 1
        
        # set voxel to start
        ix = startVoxel[0]
        iy = startVoxel[1]
        iz = startVoxel[2]
        
        #// directional constants used in scan path
        delta_ix = 1 if rayVector[0] >= 0 else -1	# determine sign
        d_x = 0 if delta_ix < 0 else 1
        delta_iy = 1 if rayVector[1] >= 0 else -1	# determine sign
        d_y = 0 if delta_iy < 0 else 1
        delta_iz = 1 if rayVector[2] >= 0 else -1	# determine sign
        d_z = 0 if delta_iz < 0 else 1
        
        #// distance from start voxel
        Dx = d_x-ix
        Dy = d_y-iy
        Dz = d_z-iz
        
        #// prevent divisions inside inner loops, so pre-calculate
        inv_e_x = 1/rayVector[0]
        inv_e_y = 1/rayVector[1]
        inv_e_z = 1/rayVector[2]
        

        #// compute line integral;
        while (ix >= 0) and (iy >= 0) and (iz >= 0) and (ix < mapDim[1]-1) and (iy < mapDim[0]-1) and (iz < mapDim[2]-1):
            
            #// possible intersections of path with voxel boundary
            #// (x, y boundary)
            #// direction of the path is taken into account
            s_x = (ix+Dx)*inv_e_x		#//s_x is total path to x-intersection
            s_y = (iy+Dy)*inv_e_y		#//s_y is total path to y-intersection
            s_z = (iz+Dz)*inv_e_z		#//s_z is total path to z-intersection
           
            # % calculate distance to each voxel boundary
            # a_x = abs((ix-x)*inv_e_x);
            # a_y = abs((iy-y)*inv_e_y);
            # a_z = abs((iz-z)*inv_e_z);
            
            
            #// only the closest intersection is really encountered
            #// find this intersection and update voxel index for this direction
            #// in some cases more boundaries are possible (45 degree paths)
            #TDOD: check influence of x voxel preference
            
            
            
            if s_x <= s_y and s_x <= s_z:				#// intersection at x-boundary
               	s_next = s_x
               	ix += delta_ix			#// x-index next voxel
                if s_x == s_y: iy += delta_iy 
                if s_x == s_z: iz += delta_iz
           
            elif s_y <= s_x and s_y <= s_z:				#// intersection at y-boundary
                s_next = s_y
                iy += delta_iy			#// y-index next voxel
                if s_y == s_z: iz += delta_iz
           
            else: #s_z <= s_x and s_z <= s_y
                s_next = s_z
                iz += delta_iz
           
           
           
            #TODO: proper length through voxel calculation
            
            
            # % calculate the length through current voxel
            # l = (((1+(e_z^2))*(a^2))^.5)*voxel_dim;
            l = np.sqrt( (1 + (rayVector[2]**2)) * (s_next**2) ) * voxel_dim
           
            # % get the attenuation coefficient
            # if AttenuationMap(round(x),round(y),round(z)) > 0
            #     tissue = AttenuationMap(round(x),round(y),round(z));
            #     atn_coef = attenuation_table(energy_index,tissue);
            # else
            #     atn_coef = 0;
            # end
            
            # LineIntegralAtt = LineIntegralAtt + l*atn_coef;
            # LineIntegralWeight = LineIntegralWeight + l*ActivityMap(round(x),round(y),round(z));

            #// update the value of the pixel  
            
            AttenuationMap[iy, ix] += sinogramVal
            
            # add value to density map, is used to normalise the output image
            DensityMap[iy, ix] += densityMapStrength
    
            #// update voxelcount and current position
            s_curr = s_next
    
        #// end of while loop for ray tracing
    
        return AttenuationMap, DensityMap