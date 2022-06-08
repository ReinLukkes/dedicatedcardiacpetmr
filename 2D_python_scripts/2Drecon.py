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

import raytrace2D

# Simulate measurement
# This is done by makeing a forward projection and then adding noise
def SimulateMeasurements(real, detector):
    measured = detector.sfp(real)
    measured = measured/np.max(measured)*100
    measured = np.random.poisson(measured)
    norm_sum = np.sum(measured)/len(detector.angles)
    
    plt.figure()
    plt.imshow(measured,cmap='gray')
    plt.show()
    
    return measured, norm_sum  

def SimulateMeasurementsC(real, detector1, detector2):
    measured = detector1.sfp(real)
    measured = measured + detector2.sfp(real)
    measured = measured/np.max(measured)*100
    measured = np.random.poisson(measured)
    norm_sum = np.sum(measured)/len(detector1.angles + detector2.angles)
    
    plt.figure()
    plt.imshow(measured,cmap='gray')
    plt.show()
    
    return measured, norm_sum  

def mlemStep(guess, measured, norm_sum, detector):
    
    # forward project guess volume
    simulated = detector.sfp(guess)
    # prevent dividing by zero
    simulated[simulated == 0] = 1e-15
    # calculate error sinogram
    error_sin = measured/simulated
    # backproject to error volume
    error_vol = detector.sbp(error_sin, N)
    # update guess volume
    guess = guess*error_vol
    # normalize guess volume
    guess = guess/np.sum(guess)*norm_sum
    
    plt.figure()
    plt.imshow(guess, cmap='gray')
    plt.show()
    
    return guess
    
def mlemStepC(guess, measured1, measured2, norm_sum, detector1, detector2):
    # forward project guess volume
    simulated1 = detector1.sfp(guess)
    simulated2 = detector2.sfp(guess)
    # prevent dividing by zero
    simulated1[simulated1 == 0] = 1e-15
    simulated2[simulated2 == 0] = 1e-15
    # calculate error sinogram
    error_sin1 = measured1/simulated1
    error_sin2 = measured2/simulated2
    # backproject to error volume
    error_vol = detector1.sbp(error_sin1, N)
    
    fig = plt.figure()
    fig.suptitle('1')
    plt.imshow(error_vol, cmap='gray')
    plt.show()
    
    error_vol += detector2.sbp(error_sin2, N)

    fig = plt.figure()
    fig.suptitle('2')
    plt.imshow(error_vol, cmap='gray')
    plt.show()
    
    # update guess volume
    guess = guess*error_vol

    fig = plt.figure()
    fig.suptitle('3')
    plt.imshow(guess, cmap='gray')
    plt.show()
    
    # normalize guess volume
    guess = guess/np.max(guess)*norm_sum
    

    fig = plt.figure()
    fig.suptitle('4')
    plt.imshow(guess, cmap='gray')
    plt.show()
    
    return guess
    
#%% When started as standalone 
if __name__ == "__main__":
    
    # Define small number to prevent dividing by zero in calculated error sinogram (see below)
    small_nr = 1e-15
    
    runSkimage      = False
    runDedicated    = True
    runFull         = True
    runCombined     = True
    
    # Import image
    real0 = imread("big_brain.png")[:,:,0]
    
    # Size of image
    N=real0.shape[0]
    guess0 = np.ones([N,N])
    guess1 = np.ones([N,N])
    guess2 = np.ones([N,N]) 
    guess3 = np.ones([N,N])  
    
    # Nr of iterations used during reconstruction
    iNrIterations = 1
    
    # define angles
    nr_angles = 10
    angles = np.arange(0, 360, 360/nr_angles)
    
    dedicated = raytrace2D.Detector(100, 100, 0, np.arange(30, 150, 120/nr_angles))
    full = raytrace2D.Detector(220, 256, 1, np.arange(0, 360, 360/nr_angles))   

    # Simulate measurement
    # This is done by makeing a forward projection and then adding noise
    
    if runSkimage:
        measured0 = radon(real0,theta=angles,circle=True)    
        measured0 = measured0/np.max(measured0)*100
        measured0 = np.random.poisson(measured0)
        norm_sum0 = np.sum(measured0)/nr_angles
        
        plt.figure()
        plt.imshow(measured0,cmap='gray')
        plt.show()
    
    if runDedicated or runCombined:
        measured1, norm_sum1 = SimulateMeasurements(real0, dedicated)
    
    if runFull or runCombined:
        measured2, norm_sum2 = SimulateMeasurements(real0, full)
    
    # if runCombined:
    #     measured3, norm_sum3 = SimulateMeasurementsC(real0, dedicated, full)
    
    
#    % MLEM loop
    for iIter in range(iNrIterations):      
    
        print(iIter+1, "/", iNrIterations)
        
        
        
        if runSkimage:
            simulated0 = radon(guess0,theta=angles,circle=True)
            simulated0[simulated0 == 0] = small_nr
            error_sin0 = measured0/simulated0
            error_vol0 = iradon(error_sin0, theta=angles, output_size=N, filter=None, interpolation='linear', circle=True)
            guess0 = guess0*error_vol0
            guess0 = guess0/np.sum(guess0)*norm_sum0
            
            plt.figure()
            plt.imshow(guess0, cmap='gray')
            plt.show()
        
        if runFull:
            guess2 = mlemStep(guess2, measured2, norm_sum2, full) 

    
        if runDedicated:
            guess1 = mlemStep(guess1, measured1, norm_sum1, dedicated)
            
        if runCombined:
            guess3 = mlemStepC(guess3, measured1, measured2, norm_sum1 + norm_sum2, dedicated, full)

    
    
    # Show reconstructed image
    
    if runSkimage:
        recon0=guess0
        
        plt.figure()
        plt.imshow(recon0,cmap='gray')
        plt.show()

    if runDedicated:
        recon1=guess1
        
        plt.figure()
        plt.imshow(recon1,cmap='gray')
        plt.show()

    if runFull:
        recon2=guess2
        
        plt.figure()
        plt.imshow(recon2,cmap='gray')
        plt.show()
        
    if runCombined:
        recon3=guess3
        
        plt.figure()
        plt.imshow(recon3,cmap='gray')
        plt.show()
        


                     