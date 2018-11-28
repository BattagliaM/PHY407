#################################################################
# Name:     AudioFilter.py                                      #
# Authors:  Michael Battaglia                                   #
# Course:   Phy407                                              #
# Inst.:    Paul Kushner                                        #
# Date:     Oct 13, 2016                                        #
# Function: Program contains routines for filtering wav audio   #
#           files using Fourier transform methods.              #
#################################################################

#essential imports
from scipy.io.wavfile import read, write
import numpy as np
import matplotlib.pyplot as plt

#function: high tone filter
def filter_high(signal, freq, cutoff):
    #Take fourier transform
    sig_F = np.fft.rfft(signal)
    #cut frequencies
    sig_Fil = np.copy(sig_F)
    sig_Fil[freq>cutoff] = 0.0
    #Returning to time domain
    sig_filtered = np.fft.irfft(sig_Fil)
    #return filtered signal and fourier transform
    return sig_filtered, sig_Fil

#function: main
if __name__ == '__main__':
    #test file to be filtered
    test_file = "GraviteaTime.wav"
    
    #read the data into two stereo channels
    sample, data = read(test_file)
    channel_0 = data[:,0]
    channel_1 = data[:,1]
    N_Points = len(data)

    #sample is sampling frequency of signal
    #related temporal sample interval
    dt = 1.0/sample #s
    #length of interval
    T = N_Points*dt
    #(angular) frequency samples
    freq = np.arange(N_Points/2+1)*2*np.pi/T
    #time samples
    t = np.arange(N_Points)*dt

    #filter frequency threshold
    cutoff = 880*2*np.pi
    #cut frequencies above threshold
    ch0_filtered, ch0_Fil = filter_high(channel_0, freq, cutoff)
    ch1_filtered, ch1_Fil = filter_high(channel_1, freq, cutoff)
    
    #write back wav file with filtered data
    data_out = np.empty(data.shape, dtype = data.dtype)
    data_out[:,0] = ch0_filtered
    data_out[:,1] = ch1_filtered
    out_file = test_file.split('.')[0]+"_Filtered."+test_file.split('.')[1]
    write(out_file, sample, data_out)

    #visualize transformations

    #Plots of original data over time
    fig, ax = plt.subplots(2,sharex=True)
    ax[0].set_title("Channel data over time")
    ax[-1].set_xlabel('Time (seconds)')
    ax[0].plot(t, channel_0)
    ax[0].set_ylabel('Channel 0 data')
    ax[1].plot(t, channel_1)
    ax[1].set_ylabel('Channel 1 data')
    fig.subplots_adjust(hspace=0)
    plt.show()
    
    #Plots of fft amplitudes of unfiltered data over frequencies
    fig, ax = plt.subplots(2,sharex=True)
    ax[0].set_title("Amplitutde of original coefficients")
    ax[-1].set_xlabel('Frequency')
    ax[0].plot(freq, abs(np.fft.rfft(channel_0)))
    ax[0].set_ylabel('Amplitude of channel 0')
    ax[1].plot(freq, abs(np.fft.rfft(channel_1)))
    ax[1].set_ylabel('Amplitude of channel 1')
    fig.subplots_adjust(hspace=0)
    plt.show()
    #Plots of fft amplitudes of filtered data over frequencies
    fig, ax = plt.subplots(2,sharex=True)
    ax[0].set_title("Amplitutde of filtered coefficients")
    ax[-1].set_xlabel('Frequency')
    ax[0].plot(freq, abs(ch0_Fil))
    ax[0].set_ylabel('Amplitude of channel 0')
    ax[1].plot(freq, abs(ch1_Fil))
    ax[1].set_ylabel('Amplitude of channel 1')
    fig.subplots_adjust(hspace=0)
    plt.show()
    
    #cut time array to between intervals
    t1 = 0.02 #s
    t2 = 0.05 #s
    index = np.logical_and(t>t1,t<t2)
    t_cut = t[index]
    ch0_cut = channel_0[index]
    ch1_cut = channel_1[index]
    ch0_fil_cut = ch0_filtered[index]
    ch1_fil_cut = ch1_filtered[index]
    #Plots of original data over time window
    fig, ax = plt.subplots(2,sharex=True)
    ax[0].set_title("Channel data over time 20ms-50ms")
    ax[-1].set_xlabel('Time (seconds)')
    ax[0].plot(t_cut, ch0_cut)
    ax[0].set_ylabel('Channel 0 data')
    ax[1].plot(t_cut, ch1_cut)
    ax[1].set_ylabel('Channel 1 data')
    fig.subplots_adjust(hspace=0)
    plt.show()
    #Plots of filtered data over time window
    fig, ax = plt.subplots(2,sharex=True)
    ax[0].set_title("Filtered channel data over time 20ms-50ms")
    ax[-1].set_xlabel('Time (seconds)')
    ax[0].plot(t_cut, ch0_fil_cut)
    ax[0].set_ylabel('Channel 0 data')
    ax[1].plot(t_cut, ch1_fil_cut)
    ax[1].set_ylabel('Channel 1 data')
    fig.subplots_adjust(hspace=0)
    plt.show()