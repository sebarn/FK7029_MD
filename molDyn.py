from numpy import*
import scipy as Sci
import scipy.linalg
from visual import*
import functions as fun
import time
from config import*
from time import clock
import datetime

#record current time
now = datetime.datetime.now()
time = now.strftime("%Y-%m-%d %H:%M")

#initialize simulation
fun.initializeFiles()
initPos = fun.initPos()
initVelo = fun.initVelo(initPos)

#initialize visualization with vpython
fun.initializeBox()
balls = []
for i in range(0,NUMBER_PARTICLES):
    balls.append(sphere(pos = initPos[i], radius = 0.05))
    
#perform calcuations for first integration step
eKin = (initVelo**2).sum()*0.5
initialForce,ePot = fun.computeForce(initPos)
eTot = eKin + ePot
currVelo = initVelo
currPos = initPos
currF = initialForce
scalingFactor = sqrt((TEMPERATURE * (NUMBER_PARTICLES) *3)/ (2*eKin) )
currVelo *= scalingFactor

energy = []

for i in range(0,NUMBER_STEPS):
    
    #carry out verlet integration
    nextPos,nextVelo, nextF,ePot = fun.integrate(currPos, currVelo, currF)
    currPos = nextPos
    currVelo = nextVelo
    currF = nextF
    
    #calculate energy
    eKin = (currVelo**2).sum()*0.5
    eTot = eKin + ePot    
    energy.append((eKin,ePot,eTot))
    
    ##rescale velocities, only used in Rahman parameters!   
    #if ((i % 100 == 0) and (i < EQUILIBRIUM_TIME) and (i>0)):
        #currVelo *= scalingFactor
        #print("rescaled!")
       
    #store calculations in files and visualize
    if (i % 1000 == 0):
        for q in range(0,NUMBER_PARTICLES):
            balls[q].pos=currPos[q]
            fun.writeEnergy(energy)
            energy = []
            
    if (i % 1000 == 0):
        print("Ekin:", eKin, "Epot:", ePot, "Etot:",eTot, "iteration:", i)
        
    if((i % 100 == 0) and (i>EQUILIBRIUM_TIME)):    
        fun.writePosVel(currPos,currVelo)
        

        
    

   


   
