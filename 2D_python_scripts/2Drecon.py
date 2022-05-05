#%% Test of 4D recon
# Casper Beijst
# 11-06-2020



#%% Packages

import numpy as np
import matplotlib.pyplot as plt
import math

from scipy.ndimage import shift,rotate

from skimage.io import imread
from skimage.transform import radon, rescale,iradon
from skimage import exposure

import radonworkshop
import raytrace2D



#%% When started as standalone 

if __name__ == "__main__":
    
    # Define small number to prevent dividing by zero in calculated error sinogram (see below)
    small_nr = 1e-15
    
    
    # Import image
    real0 = imread("big_brain.png")[:,:,0]
    
    # Size of image
    N=real0.shape[0]
    guess0 = np.ones([N,N])
    guess1 = np.ones([N,N])
    guess2 = np.ones([N,N])   
    # Nr of iterations used during reconstruction
    iNrIterations = 10
    
    # define angles
    nr_angles = 360
    angles0 = np.arange(0, 360, 360/nr_angles)
    angles1 = np.arange(30, 150, 180/nr_angles)
    
    # define scanner
    scanner_radius = 220
    circle_array = []
    
    # define cardiac detector
    detector_size = 1   #size in pixels of single detector cell
    detector_num = 180   #number of detector cells
    detector_height = 30
    detector_width = 180
    
    # ------------------------------
    
    
    
    measured2 = raytrace2D.ForwardProjection2(angles0, real0, detector_num, detector_width, scanner_radius)
    measured2 = measured2/np.max(measured2)*100
    measured2 = np.random.poisson(measured2)
    
    
    # ------------------------------
    # measured1 = raytrace2D.ForwardProjection(angles1, real0, detector_num, detector_size, detector_height)
    # measured1 = measured1/np.max(measured1)*100
    
        
    # plt.figure()
    # plt.imshow(measured1,cmap='gray')
    # plt.show()
    
    
    # measured1 = np.random.poisson(measured1)
   

    # Simulate measurement
    # This is done by makeing a forward projection and then adding noise
#    im_real = np.ones([N,N])
    # measured0 = radon(real0,theta=angles0,circle=True)    
    # measured0 = measured0/np.max(measured0)*100
    # measured0 = np.random.poisson(measured0)
    
    # measured0 = radonworkshop.radont(real0, detector_size, detector_num, detector_heigt, theta=angles, circle=True)
    
    
    # Show real imaged 
    # plt.figure()
    # plt.imshow(measured0,cmap='gray')
    # plt.show()
    
    # plt.figure()
    # plt.imshow(measured1,cmap='gray')
    # plt.show()
    
    plt.figure()
    plt.imshow(measured2,cmap='gray')
    plt.show()
    
    
    # Calculate normalization
    # norm_sum0 = np.sum(measured0)/nr_angles
    # norm_sum1 = np.sum(measured1)/nr_angles
    norm_sum2 = np.sum(measured2)/nr_angles
    

 
    
    
#    % MLEM loop
    for iIter in range(iNrIterations):      
    
        print(iIter, "/", iNrIterations)
        # forward project guess volume
        simulated0 = radon(guess0,theta=angles0,circle=True)
        
        simulated2 = raytrace2D.ForwardProjection2(angles0, guess0, detector_num, detector_width, scanner_radius)
        # simulated0 = radonworkshop.radont(guess0, detector_size, detector_num, detector_heigt,theta=angles,circle=True)
        # simulated1 = raytrace2D.ForwardProjection(angles1, guess1, detector_num, detector_size, detector_height)
        # prevent dividing by zero
        # simulated0[simulated0 == 0] = small_nr
        # simulated1[simulated1 == 0] = small_nr
        
        # calculate error sinogram
        # error_sin0 = measured0/simulated0
        # error_sin1 = measured1/simulated1
        error_sin2 = measured2/simulated2
        
        # backproject to error volume
        # error_vol0 = iradon(error_sin0, theta=angles0, output_size=N, filter=None, interpolation='linear', circle=True)
        # error_vol0 = radonworkshop.iradont(error_sin0, detector_size, detector_num, detector_heigt, theta=angles, output_size=N, filter_name=None, interpolation='linear', circle=True)
        # error_vol1 = raytrace2D.BackProjection(error_sin1, angles1, N, detector_num, detector_size, detector_height)
        error_vol2 = raytrace2D.BackProjection2(error_sin2, angles0, N, detector_num, detector_size, detector_width, scanner_radius)
        
        # update guess volume
        # guess0 = guess0*error_vol0
        # guess1 = guess1*error_vol1
        guess2 = guess2*error_vol2

        # plt.figure()
        # plt.imshow(guess0, cmap='gray')
        # plt.show()
        
        # plt.figure()
        # plt.imshow(guess1, cmap='gray')
        # plt.show()
        
        # maxpixel = np.unravel_index(np.argmax(guess0), np.shape(guess0))
        # print(maxpixel)
        
        
        # normalize guess volume <---- this is wrong for the new implementation?
        # guess0 = guess0/np.sum(guess0)*norm_sum0;
        # guess1 = guess1/np.sum(guess1)*norm_sum1;
        guess2 = guess2/np.sum(guess2)*norm_sum2;    
        
        plt.figure()
        plt.imshow(guess2, cmap='gray')
        plt.show()
        
        
        

    # recon0=guess0
    # recon1=guess1
    recon2=guess2
    
    # print(np.amax(recon0))
    
    # histogram equalisation
    # img_cdf, bin_centers = exposure.cumulative_distribution(recon0)
    # recon0 = np.interp(recon0, bin_centers, img_cdf)
    
    
    # Show reconstructed image
    # plt.figure()
    # plt.imshow(recon0,cmap='gray')
    # plt.show()

    # plt.figure()
    # plt.imshow(recon1,cmap='gray')
    # plt.show()

    plt.figure()
    plt.imshow(recon2,cmap='gray')
    plt.show()


    
                                    