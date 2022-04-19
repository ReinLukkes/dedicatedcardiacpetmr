#%% Test of 4D recon
# Casper Beijst
# 11-06-2020



#%% Packages

import numpy as np
import matplotlib.pyplot as plt

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
    
    # Size of image
    N=256
    guess0 = np.ones([N,N])
        
    # Nr of iterations used during reconstruction
    iNrIterations = 10
    
    # define angles
    nr_angles = 360
    angles = np.arange(0, 180, 180/nr_angles)
    
    # define cardiac detector
    detector_size = 1   #size in pixels of single detector cell
    detector_num = 60   #number of detector cells
    detector_height = 10

    
    
    # Import image
    real0 = imread("brain.png")[:,:,0]
    
    # ------------------------------
    measured0 = raytrace2D.ForwardProjection(angles, real0, detector_num, detector_size, detector_height)
    
   

    # Simulate measurement
    # This is done by makeing a forward projection and then adding noise
#    im_real = np.ones([N,N])
    # measured0 = radon(real0,theta=angles,circle=True)
    # measured0 = radonworkshop.radont(real0, detector_size, detector_num, detector_heigt, theta=angles, circle=True)
    measured0 = measured0/np.max(measured0)*100
    measured0 = np.random.poisson(measured0)
    
    # Show real imaged 
    plt.figure()
    plt.imshow(measured0,cmap='gray')
    plt.show()
    
    



    # Calculate normalization
    norm_sum0 = np.sum(measured0)/nr_angles

 
    
    
#    % MLEM loop
    for iIter in range(10):      
    
        # forward project guess volume
        # simulated0 = radon(guess0,theta=angles,circle=True)
        # simulated0 = radonworkshop.radont(guess0, detector_size, detector_num, detector_heigt,theta=angles,circle=True)
        simulated0 = raytrace2D.ForwardProjection(angles, guess0, detector_num, detector_size, detector_height)
        # prevent dividing by zero
        simulated0[simulated0 == 0] = small_nr
        
        # calculate error sinogram
        error_sin0 = measured0/simulated0
        
        # backproject to error volume
        # error_vol0 = iradon(error_sin0, theta=angles, output_size=N, filter=None, interpolation='linear', circle=True)
        # error_vol0 = radonworkshop.iradont(error_sin0, detector_size, detector_num, detector_heigt, theta=angles, output_size=N, filter_name=None, interpolation='linear', circle=True)
        error_vol0 = raytrace2D.BackProjection(error_sin0, angles, N, detector_num, detector_size, detector_height)
        
        # update guess volume
        guess0 = guess0*error_vol0

        
        # normalize guess volume <---- this is wrong for the new implementation?
        guess0 = guess0/np.sum(guess0)*norm_sum0;

        
        plt.figure()
        plt.imshow(guess0, cmap='gray')
        plt.show()
        

    recon0=guess0
    
    print(np.amax(recon0))
    
    # histogram equalisation
    # img_cdf, bin_centers = exposure.cumulative_distribution(recon0)
    # recon0 = np.interp(recon0, bin_centers, img_cdf)
    
    
    # Show reconstructed image
    plt.figure()
    plt.imshow(recon0,cmap='gray')
    plt.show()


    
                                    