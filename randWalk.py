#################################################################
# Name:     randWalk.py                                         #
# Authors:  Michael Battaglia                                   #
# Course:   Phy407                                              #
# Inst.:    Paul Kushner                                        #
# Date:     Nov 25, 2016                                        #
# Function: Program simulates random walk.                      #
#################################################################

#essential modules
import numpy as np
import matplotlib.pyplot as plt

#function: random walk in arbitrary dimensions
def randWalk(pos_0, lims, steps):
    #pos: position vector
    #lims: vector of dimension lengths
    #steps: number of steps to take

    #list of positions
    pos = np.zeros((steps+1,len(pos_0)))
    #initial position
    pos[0] = pos_0
    #take each step
    for i in range(steps):
        step = randstep(pos[i], lims)
        pos[i+1]=pos[i]+step
    #return path taken
    return pos

#function: random step in arbitrary dimensions
def randstep(pos_i, lims):
    #take a random step in a random direction with a random orientation
    step = np.zeros(len(pos_i))
    step[np.random.randint(0,len(pos_i))] = 1-2*np.random.randint(0,2)
    if any(pos_i+step==lims) or any(pos_i+step<0):
        #exited box, reroll step
        return randstep(pos_i, lims)
    else:
        #take step
        return step

#function: animated plot of 2D random walk
def D2plot(pos, animate=0.01):
    if animate:
        for i in range(len(pos)):
            plt.cla()
            plt.title('random walk')
            plt.plot(pos.T[0][:i+1],pos.T[1][:i+1])
            plt.xlabel('x position')
            plt.ylabel('y position')
            plt.draw()
            plt.pause(animate)
    else:
        plt.title('random walk')
        plt.plot(pos.T[0],pos.T[1])
        plt.xlabel('x position')
        plt.ylabel('y position')
    plt.show()

#function: main
if __name__ == '__main__':
    #size of box
    L = np.array([101,101])
    #central position
    central = np.array([50,50])
    #take a random walk of 1000 steps
    pos = randWalk(central, L, 5000)
    #plot path
    D2plot(pos, animate=False)
    
    
    
        
