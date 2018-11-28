#################################################################
# Name:     IdealGas.py                                         #
# Authors:  Michael Battaglia                                   #
# Course:   Phy407                                              #
# Inst.:    Paul Kushner                                        #
# Date:     December 1, 2016                                    #
# Function: Program simulates non-interacting quantum ideal gas #
#           using Monte Carlo methods. Evolves internal energy  #
#           and Boltzmann Distribution of particle energies.    #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt

#function: perform Monte Carlo simulation of idea gas
def MCIG(T, N, delta, E0=None):
    #create 2D array of quantum number initial states
    n = np.ones((N,3),int)
    #store energy at each time step
    energy = []
    #energy initial state
    if E0 is None:
        #ground state
        E = 3*N*np.pi*np.pi/2
    else:
        #initial state
        E = E0
    energy.append(E)
    #initiaize averaging window
    energy_ave = E
    error = 2*delta
    stepset = 100000
    
    #main loop
    while (error > delta):
        for k in range(stepset):
            #choose the particle
            i = np.random.randint(N)
            #choose direction
            j = np.random.randint(3)
            #choose orientation
            dn = 1 - 2*np.random.randint(2)
            #energy change associated with orientation
            dE = (dn*2*n[i,j]+1)*np.pi*np.pi/2
            #check if move is legal (quantum numbers > 0)
            if n[i,j]>1 or dn == 1:
                #check Metropolis-Hastings decision
                if np.random.random() < np.exp(-dE/T):
                    #accept step
                    n[i,j] += dn
                    E += dE
            #store energy
            energy.append(E)
        #calculate average over last window
        window_ave = np.mean(energy[-1-stepset:-1])
        #calculate error of average
        error = abs(window_ave - energy_ave)/window_ave
        #new average energy
        energy_ave = window_ave
    energy = np.array(energy)
    steps = len(energy)

    #return energy over time, and final particle states
    return n, energy, energy_ave, steps

#function: histogram of particle quantum number
def nHist(n):
    #calculate energy of each particle
    energy_n = np.square(n).sum(-1)
    #create histogram of energy count in 50 bins
    energy_freq, bin_edges = np.histogram(energy_n, 50)
    #x axis of histogram plot (n value at each bin)
    energy_vals = 0.5*(bin_edges[:-1]+bin_edges[1:])
    n_vals = energy_vals**0.5
    #calculate average quantum number
    n_ave = (energy_freq*n_vals).sum()/energy_freq.sum()

    #return histogram data
    return energy_freq, n_vals, n_ave

#function: main
if __name__ == '__main__':
    #controllable parameters
    T = 10.0 #temperature in natural units (kB*T)
    N = 1000 #number of particles
    delta = 0.01 #run until average energy obtained to 1% error

    #perform Monte Carlo sim
    n, energy, energy_ave, steps = MCIG(T, N, delta)
    print "Simulation completed to required accuracy in steps:",steps
    print energy_ave

    #plot energy against time step
    plt.title("Energy in Ideal Gas")
    plt.plot(energy)
    plt.ylabel("Energy")
    plt.xlabel("time step")
    plt.show()

    #histogram of final state
    energy_freq, n_vals, n_ave = nHist(n)
    print n_ave
    
    #plot histogram of particle energies in final state
    plt.title("Energy Distribution of Particles")
    plt.bar(n_vals, energy_freq, width=0.1)
    plt.ylabel("particles")
    plt.xlabel("quantum number")
    plt.show()

    #perform trials over a number of temperatures
    T_trials = [10.0, 40.0, 100.0, 400.0, 1200.0, 1600.0]
    E_trials = []
    n_trials = []
    #main loop
    for T in T_trials:
        #perform Monte Carlo sim
        n, energy, energy_ave, steps = MCIG(T, N, delta)
        print "Simulation completed to required accuracy in steps:",steps
        #histogram of final state
        energy_freq, n_vals, n_ave = nHist(n)

        #record trial average E and n values
        E_trials.append(energy_ave)
        n_trials.append(n_ave)
    E_trials = np.array(E_trials)
    n_trials = np.array(n_trials)

    #plot temperature relationships
    plt.title("Equilibrium Energy vs Temperature")
    plt.plot(T_trials, E_trials)
    plt.ylabel("Energy")
    plt.xlabel("Temperature")
    plt.show()

    plt.title("Average Particle Energy vs Temperature")
    plt.plot(T_trials, n_trials)
    plt.ylabel("Q number")
    plt.xlabel("Temperature")
    plt.show()
    
