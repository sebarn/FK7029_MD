import numpy as np
import scipy as sci
import functions as fun
import visual as vis
from config import*

fileNamePosition = "MC_500_positions.txt"
fileNameEnergy = "MC_500_energies.txt"
fun.initializeFilesMC(fileNamePosition,fileNameEnergy)

positions = np.zeros((NUMBER_PARTICLES,3))
positions = fun.initPos()
accepted = 0.
rejected = 0.
monitorRatio = []

currentEnergy = 0.
trialEnergy = 0.

#initialize visualization with vpython
fun.initializeBox()
balls = []
for i in range(0,NUMBER_PARTICLES):
    balls.append(vis.sphere(pos = positions[i], radius = 0.05))    

#initialize storage array for energy
energy = np.zeros(NUMBER_STEPS)

#start Metropolis Monte Carlo
for i in range(0,NUMBER_STEPS):
    trialAtom = np.random.random_integers(0,NUMBER_PARTICLES-1)
    trialMove = np.array([2*np.random.ranf() - 1.,
                          2*np.random.ranf() - 1.,
                          2*np.random.ranf() - 1.])*DJUMP
    
    currentEnergy = fun.computeEnergy(positions)
    positions[trialAtom] += trialMove
    positions[trialAtom] = fun.checkPositionsMC(positions[trialAtom])
    trialEnergy = fun.computeEnergy(positions)
    deltaEnergy = trialEnergy - currentEnergy
    
    #if energy change is negative, then accept trial configuration
    if deltaEnergy < 0.:
        accepted += 1
        #sample the trial energy
        energy[i] = trialEnergy
    #if energy change is positive, accept or reject according to metropolis
    #alogorithm (todo reference)
    elif np.random.ranf() < np.exp(-BETA * deltaEnergy):
        accepted += 1.
        #sample the trial energy
        energy[i] = trialEnergy
    else:
        #return to old configuration, if move is rejected
        positions[trialAtom] -= trialMove
        positions[trialAtom] = fun.checkPositionsMC(positions[trialAtom])
        #sample the old energy
        energy[i] = currentEnergy
        rejected += 1.
    
    #take samples of positions
    if (i % 100 == 0 and i > 10**5):
        fun.writePosMC(positions,fileNamePosition)
    
    #visualization and terminal output
    if(i% 1000 == 0 and i > 0):
        for q in range(0,NUMBER_PARTICLES):
            balls[q].pos=positions[q]
        ratio = accepted/(accepted+rejected)
        monitorRatio.append(ratio)
        accepted = 0.
        rejected = 0.
        print ratio, currentEnergy

#write energies to file    
np.savetxt(fileNameEnergy,energy,delimiter=',')
np.savetxt("ratio_MC_108.txt",monitorRatio,delimiter=',')
        