import numpy as np
from config import*
import matplotlib.pyplot as plt
import functions as func
import csv

#analyse data from Monte Carlo experiment wiht 108 particles
computeRDF = "no";

fileNameEnergy = "MC_108_energies.txt"
fileNamePosition = "MC_108_positions.txt"
fileNameRDF = "MC_108_RDF.txt"
fileNameRDFerrors = "MC_108_RDF_errors.txt"

#import energy data
energy = np.genfromtxt(fileNameEnergy,delimiter = ',', skip_header=3)

#plot energy evolution for entire time interval
plt.figure(1)
x = np.linspace(0,len(energy),len(energy))
x *= 10**(-5)
plt.ylim(-6,0)
plt.xlim(0,max(x))
plt.xlabel(r'number of steps x $10^5$')
plt.ylabel(r'$\frac{U}{N}$', fontsize=20)
plt.plot(x,energy/NUMBER_PARTICLES, label="potential")
plt.savefig('/home/sebastian/Dropbox/msc/FK7029/molDyn/tex/figures/energy_MC_108.png')

#plot energy evolution during equilibration
zoom = 10**5
plt.figure(2)
x = np.linspace(0,zoom,zoom)
x *= 10**(-5)
plt.ylim(-6,0)
plt.xlabel(r'number of steps x $10^5$')
plt.ylabel(r'$\frac{U}{N}$', fontsize=20)
plt.plot(x,energy[0:zoom]/NUMBER_PARTICLES, label="potential")
plt.savefig('/home/sebastian/Dropbox/msc/FK7029/molDyn/tex/figures/energy_MC_108_zoom.png')

#now estimate sigma for potential energy
potentialEnergy = energy[10**5::]/NUMBER_PARTICLES
sigma = func.computeSigmaSquared(potentialEnergy)
y = sigma[:,1]/sigma[:,2]
yerror = np.sqrt(2*(sigma[:,1])**2/(sigma[:,2])**3)
#save for heat capacity calculation
varPotE = y
deltaVarPotE = yerror
x = np.linspace(0,len(y)-1,len(y))
y = np.sqrt(y)
yerror = np.sqrt(yerror)
dummy = sigma

plt.figure(3)
plt.xlabel("M")
plt.xlim(0,len(x) +1)
plt.ylabel(r'$\sigma$', fontsize=16)
plt.errorbar(x,y,yerr=yerror, fmt='s', label="potential")
plt.legend(loc=2)
plt.savefig('/home/sebastian/Dropbox/msc/FK7029/molDyn/tex/figures/sigmas_MC_108.png')

#compute RDF from data, if not done previously compute it
if(computeRDF == "yes"):
    #import position data and set paramaters to calculate rad
    positions = np.genfromtxt(fileNamePosition,delimiter = ',',skip_header=3)
    #find total number of samples
    numberSamples = len(positions[:,0])/NUMBER_PARTICLES
    print("number samples:")
    print(numberSamples)
    numberBins = 300
    firstSample = numberSamples - 5000
    lastSample = numberSamples - 1000   
    
    g,u,delg,errors = func.radialDistributionFunction(positions,
                                                      numberBins,
                                                      firstSample,
                                                      lastSample)
    
    
    np.savetxt(fileNameRDF,g,delimiter=',')
    np.savetxt(fileNameRDFerrors,errors,delimiter=',')
    print(u)
else:
    #g is the radial distribution function, or g as in FS or AT
    g = np.genfromtxt(fileNameRDF,delimiter=',');
    errors = np.genfromtxt(fileNameRDFerrors, delimiter=',');
    delg = HALF_BOX_LENGTH/float(len(g));

#select points for error analysis in RDF 126 (first peak), 182 (first mininum)
#243 (second peak), 295(tail)
#plot RDF

#normalize errors
nhis = 300
for i in range(0,nhis):
		vb = (float((i+1)**3) - float(i**3)) * delg**3
		nid = 4./3. * np.pi * vb * DENSITY
		errors[:,i] = errors[:,i] / (float(NUMBER_PARTICLES)*nid)
#select four interesting points from errors		
errors = errors[:,(126,182,243,295)]
errors = np.concatenate((np.zeros((1,4)),errors),axis=0)
errors = [np.diff(errors[:,k]) for k in range(4)]
#use Flyvberg&Peterson method to estimate sigma
errors = [func.computeSigmaSquared(errors[k]) for k in range(4)]

plt.figure(4)
plotTitles = ["first peak","first minimum","second peak", "tail"];
j = 1
for j in range(0,4):
	sigma = errors[j]
	y = sigma[:,1]/sigma[:,2]
	yerror = np.sqrt(2*(sigma[:,1])**2/(sigma[:,2])**3)
	x = np.linspace(0,len(y)-1,len(y))
	y = np.sqrt(y)
	yerror = np.sqrt(yerror)
	
	plt.subplot(2,2,j+1)
	plt.title(plotTitles[j])
	plt.xlabel("M")
	plt.xlim(0,len(x))
	plt.ylabel(r'$\sigma$', fontsize=16)
	plt.errorbar(x,y,yerr=yerror, fmt='s')
#plt.savefig('/home/sebastian/Dropbox/msc/FK7029/molDyn/tex/figures/errors_RDF_MC_108.png')

#this I have to do by hand, because I cannot automatize
#reading the correct sigma from the plots
yerror = np.zeros(nhis)
xError = np.array([126,182,243,295])*delg
yError = np.array([g[126],g[182],g[243],g[295]])
yErrorBar = np.array([0.0237,0.0078,0.0010, 0.0062])
	
#plot the RDF and the errors just calculated
plt.figure(5)
plt.xlabel('r')
plt.ylabel('g(r)')
x = np.linspace(0,len(g)-1,len(g))*delg
plt.plot(x,g)
plt.errorbar(xError,yError,yerr=yErrorBar, fmt='x')
plt.savefig('/home/sebastian/Dropbox/msc/FK7029/molDyn/tex/figures/RDF_MC_108.png')

def heatCapacity(TEMPERATURE,var):
 result = 1.5 + var/TEMPERATURE**2
 return result

#calculate heat capacity and error
#temperature
ePot = energy[10000::]
#estimate the sample length roughly from Flyvberg&Peterson method
sampleLength = 200000
heatCap = np.zeros(10)

k =0
for  n in range(0,10):
 sample = ePot[n*sampleLength:(n+1)*sampleLength]
 varSample = np.var(sample)/NUMBER_PARTICLES
 heatCap[n] = heatCapacity(TEMPERATURE,varSample)
