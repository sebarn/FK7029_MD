import numpy as np
import scipy as Sci
import scipy.linalg
from config import*
from visual import*
from scipy.stats import maxwell
from scipy import weave
from scipy.weave import converters
from time import clock
import datetime

def initPos():
    b = LATTICE_CONSTANT
    x0 = ones((4,3))
    x0 *= b / 2.
    x0[0,:] = 0.
    x0[1,2] = 0.
    x0[2,1] = 0.
    x0[3,0] = 0.
    positions = np.zeros((NUMBER_PARTICLES,3))    
    index = 0    
    for i in range(L):
        for j in range(L):
            for k in range(L):
                for chi in range(4):
                    positions[index,0] = x0[chi,0] + i*b
                    positions[index,1] = x0[chi,1] + j*b
                    positions[index,2] = x0[chi,2] + k*b
                    index += 1    
    positions += b*0.5
    return positions


def initVelo(currentPositions):
    #get numper of particles from main routine and create velocities accordingly
    
    totVelocity = array([0.,0.,0.])
    kineticEnergy = 0.
    particleVelocities = np.zeros((NUMBER_PARTICLES,3))
    #draw initial velocities from normal distribution
    mu = 0.
    sigma = 0.5
    
    for component in range(0,3):
        particleVelocities[:,component] = scipy.stats.maxwell.rvs(size=NUMBER_PARTICLES)
        averageVelocity = particleVelocities[:,component].sum()/NUMBER_PARTICLES
        particleVelocities[:,component] = particleVelocities[:,component] - averageVelocity
      
    totVelocity = particleVelocities.sum()
            
    print("total veolcity (should be zero):", totVelocity)
    return particleVelocities



def integrate(currPos, currVel, currF):
    nextPos = currPos + TIME_STEP * currVel + TIME_STEP*TIME_STEP*0.5*currF
    nextPos = checkParticlePosition(nextPos)
    nextF,ePot = computeForce(nextPos)
    nextVel = currVel + TIME_STEP * 0.5 *(currF + nextF)
    return (nextPos, nextVel, nextF, ePot)    

def checkParticlePosition(positions):
    code = """
    using namespace blitz;
        
    for(int k = 0; k < NUMBER_PARTICLES; k++){
        for(int j = 0; j < 3; j ++){
            if(positions(k,j) > BOX_LENGTH){
                positions(k,j) -= BOX_LENGTH;
                }
            else if(positions(k,j) < 0.){
                positions(k,j) += BOX_LENGTH;
                }
            }//end inner for
        }//end outer for
    """
    
    arguments = ['positions','NUMBER_PARTICLES','BOX_LENGTH']
    weave.inline(code,arguments,type_converters = converters.blitz,compiler = 'gcc')
    return positions
def checkPositionsMC(pos):
    #maybe speed this up by using stacks
    for k in range(0,3):
            if pos[k] > BOX_LENGTH:
                pos[k] -= BOX_LENGTH
            elif pos[k] < 0.:
                pos[k] += BOX_LENGTH
    return pos

def initializeBox():
    #todo
    face1 = box(pos=(BOX_LENGTH/2,BOX_LENGTH/2,0), 
                size = (BOX_LENGTH,BOX_LENGTH,0.01),
                color=color.green, opacity = 0.2)
    face2 = box(pos=(0,BOX_LENGTH/2,BOX_LENGTH/2), 
                size = (0.01,BOX_LENGTH,BOX_LENGTH),
                color=color.green, opacity = 0.2)    
    face3 = box(pos=(BOX_LENGTH/2,0,BOX_LENGTH/2),
                size = (BOX_LENGTH,0.01,BOX_LENGTH),
                color=color.green, opacity = 0.2)    
    face4 = box(pos=(BOX_LENGTH/2,BOX_LENGTH/2,BOX_LENGTH),
                size = (BOX_LENGTH,BOX_LENGTH,0.01),
                color=color.green, opacity = 0.2)
    face5 = box(pos=(BOX_LENGTH,BOX_LENGTH/2,BOX_LENGTH/2),
                size = (0.01,BOX_LENGTH,BOX_LENGTH),
                color=color.green, opacity = 0.2)    
    face6 = box(pos=(BOX_LENGTH/2,BOX_LENGTH,BOX_LENGTH/2), 
                size = (BOX_LENGTH,0.01,BOX_LENGTH),
                color=color.green, opacity = 0.2)        



#ffunctions for force and energy calculations 
def computeForce(currentPositions):        
    code = """
    using namespace blitz;
    Array<double,1> distance(3);
    double distanceSquared, r2i, r6i, lennardJones;
    double potentialEnergy = 0.;
    
    for( int iParticle = 0; iParticle < (NUMBER_PARTICLES - 1); iParticle++){
        for( int jParticle = iParticle + 1; jParticle < NUMBER_PARTICLES; jParticle++){
            distance(0) = currentPositions(iParticle,0)-currentPositions(jParticle,0);
            distance(0) = distance(0) - BOX_LENGTH * round(distance(0)/BOX_LENGTH);
            distance(1) = currentPositions(iParticle,1)-currentPositions(jParticle,1);
            distance(1) = distance(1) - BOX_LENGTH * round(distance(1)/BOX_LENGTH);
            distance(2) = currentPositions(iParticle,2)-currentPositions(jParticle,2);
            distance(2) = distance(2) - BOX_LENGTH * round(distance(2)/BOX_LENGTH);
            distanceSquared = distance(0)*distance(0) + distance(1)*distance(1) + distance(2)*distance(2);
            if(distanceSquared < CUT_OFF_RADIUS_SQUARED){
                r2i = 1./distanceSquared;
                r6i = r2i * r2i * r2i;
                lennardJones = 48. * r2i * r6i * (r6i - 0.5);
                force(iParticle,0) += lennardJones*distance(0);
                force(iParticle,1) += lennardJones*distance(1);
                force(iParticle,2) += lennardJones*distance(2);
                force(jParticle,0) -= lennardJones*distance(0);
                force(jParticle,1) -= lennardJones*distance(1);
                force(jParticle,2) -= lennardJones*distance(2);
                potentialEnergy += 4.* r6i * (r6i - 1.) - CUT_OFF_ENERGY;
                                
                }
            
            }//end inner for loop
    }//end outer for loop
    return_val = potentialEnergy;
    
    """
    #args that are passed into weave.inline and created inside computeForce
    #potentialEnergy = 0.
    force = np.zeros((NUMBER_PARTICLES,3))
    
    #all args
    arguments = ['currentPositions','force','NUMBER_PARTICLES','CUT_OFF_RADIUS_SQUARED','BOX_LENGTH','CUT_OFF_ENERGY']
    #evaluate stuff in code
    potentialEnergy = weave.inline(code,arguments,type_converters = converters.blitz,compiler = 'gcc')    
      
    return force, potentialEnergy

def computeEnergy(currentPositions):        
    code = """
    using namespace blitz;
    Array<double,1> distance(3);
    double distanceSquared, r2i, r6i, lennardJones;
    double potentialEnergy = 0.;
    
    for( int iParticle = 0; iParticle < (NUMBER_PARTICLES - 1); iParticle++){
        for( int jParticle = iParticle + 1; jParticle < NUMBER_PARTICLES; jParticle++){
            distance(0) = currentPositions(iParticle,0)-currentPositions(jParticle,0);
            distance(0) = distance(0) - BOX_LENGTH * round(distance(0)/BOX_LENGTH);
            distance(1) = currentPositions(iParticle,1)-currentPositions(jParticle,1);
            distance(1) = distance(1) - BOX_LENGTH * round(distance(1)/BOX_LENGTH);
            distance(2) = currentPositions(iParticle,2)-currentPositions(jParticle,2);
            distance(2) = distance(2) - BOX_LENGTH * round(distance(2)/BOX_LENGTH);
            distanceSquared = distance(0)*distance(0) + distance(1)*distance(1) + distance(2)*distance(2);
            if(distanceSquared < CUT_OFF_RADIUS_SQUARED){
                r2i = 1./distanceSquared;
                r6i = r2i * r2i * r2i;
                lennardJones = 48. * r2i * r6i * (r6i - 0.5);
                potentialEnergy += 4.* r6i * (r6i - 1.) - CUT_OFF_ENERGY;
                                
                }
            
            }//end inner for loop
    }//end outer for loop
    return_val = potentialEnergy;
    
    """
    #args that are passed into weave.inline and created inside computeForce
    #potentialEnergy = 0.
    force = np.zeros((NUMBER_PARTICLES,3))
    
    #all args
    arguments = ['currentPositions','NUMBER_PARTICLES','CUT_OFF_RADIUS_SQUARED','BOX_LENGTH','CUT_OFF_ENERGY']
    #evaluate stuff in code
    potentialEnergy = weave.inline(code,arguments,type_converters = converters.blitz,compiler = 'gcc')    
      
    return potentialEnergy



#functions for file input and output
def initializeFiles():
    now = datetime.datetime.now()
    time = now.strftime("%Y-%m-%d %H:%M")
    #initialize file for positions data
    FILE = open("positions.txt", "a")
    FILE.write(time)
    FILE.write("\n")
    FILE.write("NUMBER_PARTICLES: ")
    FILE.write(str(NUMBER_PARTICLES))
    FILE.write("\t TIME_STEP: ")
    FILE.write(str(TIME_STEP))
    FILE.write("\t NUMBER_STEPS: ")
    FILE.write(str(NUMBER_STEPS))
    FILE.write("\t TEMPERATURE: ")
    FILE.write(str(TEMPERATURE))
    FILE.write("\n \n")
    FILE.close()
    #initialize file for velocity data
    FILE = open("velocities.txt", "a")
    FILE.write(time)
    FILE.write("\n")
    FILE.write("NUMBER_PARTICLES: ")
    FILE.write(str(NUMBER_PARTICLES))
    FILE.write("\t TIME_STEP: ")
    FILE.write(str(TIME_STEP))
    FILE.write("\t NUMBER_STEPS: ")
    FILE.write(str(NUMBER_STEPS))
    FILE.write("\t TEMPERATURE: ")
    FILE.write(str(TEMPERATURE))
    FILE.write("\n \n")
    FILE.close()
    #same for energy
    FILE = open("energy.txt", "a")
    FILE.write(time)
    FILE.write("\n")
    FILE.write("NUMBER_PARTICLES: ")
    FILE.write(str(NUMBER_PARTICLES))
    FILE.write("\t TIME_STEP: ")
    FILE.write(str(TIME_STEP))
    FILE.write("\t NUMBER_STEPS: ")
    FILE.write(str(NUMBER_STEPS))
    FILE.write("\t TEMPERATURE: ")
    FILE.write(str(TEMPERATURE))
    FILE.write("\n \n")
    FILE.close()    

def writePosVel(positions, velocities):
    FILE = open("positions.txt","a")
    np.savetxt(FILE,positions,delimiter=',')    
    FILE.close()
    FILE = open("velocities.txt","a")
    np.savetxt(FILE,velocities,delimiter=',')
    FILE.close()
def writeEnergy(energy):
    FILE = open("energy.txt","a")
    np.savetxt(FILE,energy,delimiter=',')
    FILE.close()
    
def initializeFilesMC(fileNameEnergy,fileNamePosition):
    now = datetime.datetime.now()
    time = now.strftime("%Y-%m-%d %H:%M")
    
    #initialize file for positions data
    FILE = open(fileNamePosition, "a")
    FILE.write(time)
    FILE.write("\n")
    FILE.write("NUMBER_PARTICLES: ")
    FILE.write(str(NUMBER_PARTICLES))
    FILE.write("\t TIME_STEP: ")
    FILE.write(str(TIME_STEP))
    FILE.write("\t NUMBER_STEPS: ")
    FILE.write(str(NUMBER_STEPS))
    FILE.write("\t TEMPERATURE: ")
    FILE.write(str(TEMPERATURE))
    FILE.write("\n \n")
    FILE.close()
    
    #initialize file for energy data
    FILE = open(fileNameEnergy, "a")
    FILE.write(time)
    FILE.write("\n")
    FILE.write("NUMBER_PARTICLES: ")
    FILE.write(str(NUMBER_PARTICLES))
    FILE.write("\t TIME_STEP: ")
    FILE.write(str(TIME_STEP))
    FILE.write("\t NUMBER_STEPS: ")
    FILE.write(str(NUMBER_STEPS))
    FILE.write("\t TEMPERATURE: ")
    FILE.write(str(TEMPERATURE))
    FILE.write("\n \n")
    FILE.close()
    
def writePosMC(positions,fileNamePosition):
    FILE = open(fileNamePosition,"a")
    np.savetxt(FILE,positions,delimiter=',')    
    FILE.close()
    
#functions used for data analysis start here
def computeSigmaSquared(qtty):
    #computes an estimate of the variance of input quantitity using method by 
    #Flyvberg & Peterson described in Frenkel/Smit Appendix D.3 input is the
    #quantitiy, for instance potential energy as numpy array
    dataIterator = 0
    data = np.zeros((100,3))
    data[dataIterator,0] = np.mean(qtty)
    data[dataIterator,1] = np.var(qtty)
    data[dataIterator,2] = len(qtty)
    dataIterator = 1
    
    #start the boxing
    length = len(qtty)
    length = floor_divide(length,2)
    box = np.zeros(length)
    while len(box) > 2:
        qttyIterator = 0
        boxIterator = 0
        while boxIterator < length:
            box[boxIterator] = 0.5 * (qtty[qttyIterator] + qtty[qttyIterator+1])
            qttyIterator += 2
            boxIterator += 1
        data[dataIterator,0] = np.mean(box)
        data[dataIterator,1] = np.var(box)
        data[dataIterator,2] = len(box)
        dataIterator += 1
        #qtty becomes old box
        qtty = box
        #cut box into half
        length = np.floor_divide(length,2)
        box = np.zeros(length)
    #select nonzero elements and return them
    data = data[0:dataIterator,:]
    return data

def potentialLJ(distance):
	value = 4. * ((1/distance)**12 - (1/distance)**6) - CUT_OFF_ENERGY
	return value

#radial distribution function
def radialDistributionFunction(data,nhis,a,b):
	#calculates the radial distribution function (RDF). input nhis is the number
	#of histograms. a is the lower limit of the range of samples, b is the upper
	#limit of the range of samples. the code is based on Frenkel/Smit, p.86
	ngr = 0
	delg = BOX_LENGTH / (2*nhis)
	g = np.zeros(nhis)
	N = NUMBER_PARTICLES
	
	#this is for the error calculation
	errors = np.zeros((b-a,nhis),dtype=np.int)
	#use this counter to fill up errors in for loop below
	errorCounter = 0
	
	for n in range(a,b):	
		pos = data[n*N:(n+1)*N]
		ngr = ngr + 1
		
		
		for i in range(0, N - 1):
			for j in range(i+1, N):
				xr = pos[i] - pos[j]
				xr = xr - BOX_LENGTH * (xr/BOX_LENGTH).round()
				r = np.linalg.norm(xr)
				if r < HALF_BOX_LENGTH:
					ig = int(r/delg)
					g[ig] += 2
		errors[errorCounter,:] = g
		errorCounter += 1
	
	for i in range(0,nhis):
		r = delg*float(i+0.5)
		vb = (float((i+1)**3) - float(i**3)) * delg**3
		nid = 4./3. * np.pi * vb * DENSITY
		#normalize g, i.e. include number of iterations
		g[i] = g[i] / (float(ngr)*float(N)*nid)
		#normalize error, i.e. exclude number of iterations
		#errors[:,i] = errors[:,i] / (float(N)*float(nid))
	
	u = 0
	for bin in range(0,nhis):
		dist = delg*(bin + 0.5)
		u += dist**2 * potentialLJ(dist) * g[bin] * delg
	
	u *= 2*np.pi*DENSITY	
	
	return g,u,delg,errors