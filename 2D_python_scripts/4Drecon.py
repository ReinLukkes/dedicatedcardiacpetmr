#%% Test of 4D recon
# Casper Beijst
# 11-06-2020



#%% Packages

import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import shift,rotate

from skimage.io import imread
from skimage.transform import radon, rescale,iradon


#%% When started as standalone 

if __name__ == "__main__":
    
    plt.close('all')
    small_nr = 1e-15
    
    N=256
    guess0 = np.ones([N,N])
        
    
    iNrIterations = 2
    
    nr_angles = 120.
    angles = np.arange(0,360,360/nr_angles)
    
    
    real0 = imread("brain.png")[:,:,0]
    real1 = rotate(real0,-10,reshape =False)
    real1 = shift(real1,[0,10],mode='nearest')
    
    real2 = rotate(real0,-20,reshape =False)
    real2 = shift(real2,[0,20],mode='nearest')
    
#    im_real = np.ones([N,N])
    measured0 =  radon(real0,theta=angles,circle=True)
    measured0 = measured0/np.max(measured0)*100
    measured0 = np.random.poisson(measured0)
    
    measured1 =  radon(real1,theta=angles,circle=True)
    measured1 = measured1/np.max(measured1)*100
    measured1 = np.random.poisson(measured1)
    
    measured2 =  radon(real2,theta=angles,circle=True)
    measured2 = measured2/np.max(measured2)*100
    measured2 = np.random.poisson(measured2)
    
    plt.figure()
    plt.imshow(real0,cmap='gray')
    plt.show()
    
    plt.figure()
    plt.imshow(real1,cmap='gray')
    plt.show()
    
    plt.figure()
    plt.imshow(real2,cmap='gray')
    plt.show()
    
#    plt.figure()
#    plt.imshow(measured0)
#    plt.show()
#    
#    plt.figure()
#    plt.imshow(measured1)
#    plt.show()
#    
    measured = measured0+measured1+measured2

    norm_sum0 = np.sum(measured0)/nr_angles
    norm_sum1 = np.sum(measured1)/nr_angles
    norm_sum2 = np.sum(measured2)/nr_angles
 
    
    
#    % MLEM
    for iIter in range(iNrIterations):
#        for i in nrtimessets        
    
        # forward project guess volume
        print(iIter)
        simulated0 = radon(guess0,theta=angles,circle=True)
        
        
        # prevent dividing by zero
        simulated0[simulated0 == 0] = small_nr
        
        # calculate error sinogram
        error_sin0 = measured0/simulated0
        
        # backproject to error volume
        error_vol0 = iradon(error_sin0, theta=angles, output_size=N, filter=None, interpolation='linear', circle=True)
        
        # apply joint objective function?
        
        
        # update guess volume
        guess0 = guess0*error_vol0
        
        # normalize guess volume
        guess0 = guess0/np.sum(guess0)*norm_sum0;
        
        
        # shift guess
        guess1 = rotate(guess0,-10,reshape =False)
        guess1 = shift(guess1,[0,10],mode='nearest')
        
        # forward project guess volume
        simulated1 = radon(guess1,theta=angles,circle=True)
        
        
        # prevent dividing by zero
        simulated1[simulated1 == 0] = small_nr
        
        # calculate error sinogram
        error_sin1 = measured1/simulated1
        
        # backproject to error volume
        error_vol1 = iradon(error_sin1, theta=angles, output_size=N, filter=None, interpolation='linear', circle=True)
        
        # apply joint objective function?
        
        
        # update guess volume
        guess1 = guess1*error_vol1
        
        # normalize guess volume
        guess1 = guess1/np.sum(guess1)*norm_sum1;
        
        
        # Shift guess back
        guess0 = shift(guess1,[0,-10],mode='nearest')
        guess0 = rotate(guess0,10,reshape=False)
        
        # shift guess
        guess2 = rotate(guess0,-20,reshape =False)
        guess2 = shift(guess2,[0,20],mode='nearest')
        
        # forward project guess volume
        simulated2 = radon(guess2,theta=angles,circle=True)
        
        
        # prevent dividing by zero
        simulated2[simulated2 == 0] = small_nr
        
        # calculate error sinogram
        error_sin2 = measured2/simulated2
        
        # backproject to error volume
        error_vol2 = iradon(error_sin2, theta=angles, output_size=N, filter=None, interpolation='linear', circle=True)
        
        # apply joint objective function?
        
        
        # update guess volume
        guess2 = guess2*error_vol2
        
        # normalize guess volume
        guess2 = guess2/np.sum(guess2)*norm_sum2;
        
        
        # Shift guess back
        guess0 = shift(guess2,[0,-20],mode='nearest')
        guess0 = rotate(guess0,20,reshape=False)
        

    recon0=guess0
    
    recon1=  rotate(recon0,-10,reshape =False)
    recon1 = shift(recon1,[0,10],mode='nearest')
    
    recon2=  rotate(recon0,-20,reshape =False)
    recon2 = shift(recon2,[0,20],mode='nearest')   
     
    plt.figure()
    plt.imshow(recon0,cmap='gray')
    plt.show()

    plt.figure()
    plt.imshow(recon1,cmap='gray')
    plt.show()
    
    plt.figure()
    plt.imshow(recon2,cmap='gray')
    plt.show()
    
