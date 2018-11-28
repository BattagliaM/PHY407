#################################################################
# Name:     integrator.py                                       #
# Authors:  Michael Battaglia                                   #
# Course:   Phy407                                              #
# Inst.:    Paul Kushner                                        #
# Date:     Sept 29, 2016                                       #
# Function: Program integrates a function numerically using     #
#           trapezoidal method, simpson's rule, and gaussian    #
#           quadrature, comparing accuracy of the three.        #
#################################################################

#essential imports
import numpy as np
import matplotlib.pyplot as plt

#function: test function to be integrated
def test_func(x):
    return 4/(1+np.square(x))

#function: integrates using trapezoidal rule
def intTrap(func, x1, x2, N):
    #make partition of integration axis into N segments
    x, dx = np.linspace(x1,x2,N+1,endpoint=True,retstep=True)
    #integrate func on domain using trapezoidal rule
    integral = 0.0
    for i in range(N): #loop over the N segments
        #calculate trapezoid area for each segment
        integral += 0.5*dx*(func(x[i+1])+func(x[i]))
    #return trapezoidal sum
    return integral

#function: integrates using simpson's rule
def intSimp(func, x1, x2, N):
    #take only even partitions
    N = N - N%2
    #make partition of integration axis into N subintervals
    x, dx = np.linspace(x1,x2,N+1,endpoint=True,retstep=True)
    #integrate func on domain using trapezoidal rule
    integral = 0.0
    for i in range(1,N/2+1): #take every 2 consecutive intervals
        #integrate parabola defined by 3 points in the two intervals
        integral += (dx/3)*(func(x[2*i-2])+4*func(x[2*i-1])+func(x[2*i]))
    #return trapezoidal sum
    return integral

#function: calculates gaussian quadrature between -1,1
def gaussxw(N):
    #initial approximation to roots of the Legendre polynomial
    a = np.linspace(3,4*N-1,N)/(4*N+2)
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))

    #Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = np.ones(N,float)
        p1 = np.copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))
    #calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)
    
    #return optimal samples with their respective weights
    return x,w

#function: calculates gaussian quadrature between a,b
def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

#function: integrates using gaussian quadrature
def intGaus(func, x1, x2, N):
    x,w = gaussxwab(N,x1,x2)
    return sum(w*func(x))

#function: main
if __name__ == '__main__':
    #computational parameters
    N_trials = range(10,500,10) #integration resolutions to be tried
    #integration bounds
    x1 = 0
    x2 = 1.0
    #analytic integral
    analytic = np.pi

    trap_err = np.zeros(len(N_trials),float)
    simp_err = np.zeros(len(N_trials),float)
    gaus_err = np.zeros(len(N_trials),float)
    gaus_est = np.zeros(len(N_trials),float)
    #integrate func using different methods for each trial N
    for i, N in enumerate(N_trials):
        #compute integrals
        trap_int = intTrap(test_func, x1, x2, N)
        simp_int = intSimp(test_func, x1, x2, N)
        gaus_int = intGaus(test_func, x1, x2, N)
        gaus_int2 = intGaus(test_func, x1, x2, 2*N)
        #compute error from analytic value
        trap_err[i] = abs(trap_int-analytic)/analytic
        simp_err[i] = abs(simp_int-analytic)/analytic
        gaus_err[i] = abs(gaus_int-analytic)/analytic
        gaus_est[i] = abs(gaus_int2-gaus_int)/gaus_int2
        
    #plot error against integration resolution
    plt.title("Error of Integration Methods")
    plt.ylabel("Fractional Error")
    plt.xlabel("Resolution")
    plt.plot(N_trials, trap_err, label="Trapezoidal Method")
    plt.plot(N_trials, simp_err, label="Simpson's Rule")
    plt.plot(N_trials, gaus_err, label="Gaussian Quadrature")
    plt.plot(N_trials, gaus_est, label="Gauss Error Estimate")
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.show()
