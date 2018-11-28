#################################################################
# Name:     Deconvolution.py                                    #
# Authors:  Michael Battaglia                                   #
# Course:   Phy407                                              #
# Inst.:    Paul Kushner                                        #
# Date:     Oct 13, 2016                                        #
# Function: Program contains routines for image processing      #
#           using Fourier transform methods.                    #
#################################################################

#essential imports
import numpy as np
import matplotlib.pyplot as plt

#function: generate gaussian PSF
def gauss_PSF(image, sigma):
    #create point spread function
    PSF = np.zeros(image.shape)
    for i in range(K):
        ip=i
        if ip>K/2:
            ip -= K #bottom half of rows moved to negative values
        for j in range(L):
            jp=j
            if jp>L/2:
                jp -= L #right half of columns moved to negative values
            PSF[i,j]=np.exp(-(ip**2+jp**2)/(2.0*sigma**2)) #compute gaussian
    return PSF

#function: deconvolution
def deconv(image, PSF):
    #fourier transform both image and PSF
    F_image = np.fft.rfft2(image)
    F_PSF = np.fft.rfft2(PSF)
    
    #pad fourier transform of PSF
    F_PSF[F_PSF<1.0e-3] = 1.0
    #deconvolve image by PSF
    D_image = F_image/(F_PSF)
    
    #inverse fourier transform deconvolved image
    d_image = np.fft.irfft2(D_image)
    return d_image

#function: main
if __name__ == '__main__':
    #load sample text containing image data
    sample_data = "blur.txt"
    image = np.loadtxt(sample_data)
    #image data
    L = image.shape[0]
    K = image.shape[1]
    sigma = 25 #gaussian point spread
    #display unprocessed image
    plt.imshow(image, cmap='gray')
    plt.show()
    
    #create point spread function
    PSF = gauss_PSF(image, sigma)
    #plot PSF
    #plt.imshow(PSF, cmap='gray')
    #plt.show()
    
    #deconvolve image by PSF
    d_image = deconv(image, PSF)
    #display processed image
    plt.imshow(d_image, cmap='gray')
    plt.show()