import numpy as np
from config import*
import matplotlib.pyplot as plt
import functions as func
import csv

#analyse data from Molecular Dynamics Carlo experiment wiht 108 particles
computeRDF = "no";

fileNameEnergy = "MD_500_energies.txt"
fileNamePosition = "MD_500_positions.txt"
fileNameRDF = "MD_500_RDF.txt"
fileNameRDFerrors = "MD_500_RDF_errors.txt"

#import energy data
energy = np.genfromtxt(fileNameEnergy,delimiter = ',',skip_header=3)
		

#plot energy evolution for entire time interval
plt.figure(1)
x = np.linspace(0,len(energy)-1,len(energy))
x *= 10**(-5)
plt.xlabel(r'number of steps x $10^5$')
plt.ylabel(r'$\frac{U}{N}$', fontsize=20)
plt.plot(x,energy[:,0]/NUMBER_PARTICLES, label="potential")
plt.plot(x,energy[:,1]/NUMBER_PARTICLES, label="kinetic")
plt.plot(x,energy[:,2]/NUMBER_PARTICLES, label="total")
plt.legend(loc=5)
plt.savefig('/home/sebastian/Dropbox/msc/FK7029/molDyn/tex/figures/energy_MD_500.png')

#plot energy evolution during equilibration
zoom = 5000
plt.figure(2)
x = np.linspace(0,zoom-1,zoom)
#time *= 10**(-5)
plt.xlabel('number of steps')
plt.ylabel(r'$\frac{U}{N}$', fontsize=20)
plt.plot(x,energy[0:zoom,0]/NUMBER_PARTICLES, label="kinetic")
plt.plot(x,energy[0:zoom,1]/NUMBER_PARTICLES, label="potential")
plt.plot(x,energy[0:zoom,2]/NUMBER_PARTICLES, label="total")
plt.legend(loc=5)
plt.savefig('/home/sebastian/Dropbox/msc/FK7029/molDyn/tex/figures/energy_MD_500_zoom.png')

#now estimate sigma for kinetic energy
eKin = energy[10**4::,0]/NUMBER_PARTICLES
sigma = func.computeSigmaSquared(eKin)
y = sigma[:,1]/sigma[:,2]
yerror = np.sqrt(2*(sigma[:,1])**2/(sigma[:,2])**3)
#save this for the heat capacity calculation below
varKinE = y
deltaVarKinE = yerror
x = np.linspace(0,len(y)-1,len(y))
y = np.sqrt(y)
yerror = np.sqrt(yerror)


plt.figure(3)
plt.xlabel("M")
plt.xlim(0,len(x) +1)
plt.ylabel(r'$\sigma$', fontsize=16)
plt.errorbar(x,y,yerr=yerror, fmt='s', label="kinetic")

#now the same for potential energy
ePot = energy[10**4::,1]/NUMBER_PARTICLES
sigma = func.computeSigmaSquared(ePot)
y = sigma[:,1]/sigma[:,2]
yerror = np.sqrt(2*(sigma[:,1])**2/(sigma[:,2])**3)
x = np.linspace(0,len(y)-1,len(y))
y = np.sqrt(y)
yerror = np.sqrt(yerror)
plt.errorbar(x,y,yerr=yerror, fmt='*', label="potential")
plt.legend(loc=2)
plt.savefig('/home/sebastian/Dropbox/msc/FK7029/molDyn/tex/figures/sigmas_MD_500.png')

#compute RDF from data, if not done previously compute it
if(computeRDF == "yes"):
    #import position data and set paramaters to calculate rad
    positions = np.genfromtxt(fileNamePosition,delimiter = ',',skip_header=3)
    #find total number of samples
    numberSamples = len(positions[:,0])/NUMBER_PARTICLES
    print("number samples:")
    print(numberSamples)
    numberBins = 300
    firstSample = numberSamples - 1010
    lastSample = numberSamples - 10
    
    
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

#select points for error analysis

#normalize errors
nhis = 300
for i in range(0,nhis):
		vb = (float((i+1)**3) - float(i**3)) * delg**3
		nid = 4./3. * np.pi * vb * DENSITY
		errors[:,i] = errors[:,i] / (float(NUMBER_PARTICLES)*nid)
#select four interesting points from errors		
errors = errors[:,(80,146,223,290)]
errors = np.concatenate((np.zeros((1,4)),errors),axis=0)
errors = [np.diff(errors[:,k]) for k in range(4)]
#use Flyvberg&Peterson method to estimate sigma
errors = [func.computeSigmaSquared(errors[k]) for k in range(4)]

plt.figure(4)
plotTitles = ["first peak","second peak","third peak", "tail"];
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
#plt.savefig('/home/sebastian/Dropbox/msc/FK7029/molDyn/tex/figures/errors_RDF_MD_500.png')

#this I have to do by hand, because I cannot automatize
#reading the correct sigma from the plots
yerror = np.zeros(nhis)
xError = np.array([80,146,223,290])*delg
yError = np.array([g[80],g[146],g[223],g[290]])
yErrorBar = np.array([0.008,0.003,0.002, 0.0015])
	
#plot the RDF and the errors just calculated
plt.figure(5)
plt.xlabel('r')
plt.ylabel('g(r)')
x = np.linspace(0,len(g)-1,len(g))*delg
plt.plot(x,g)
plt.errorbar(xError,yError,yerr=yErrorBar, fmt='x')
plt.savefig('/home/sebastian/Dropbox/msc/FK7029/molDyn/tex/figures/RDF_MD_500.png')

def heatCapacity(T,var):
 result = 1./( (2./3.) - ((4.* var) / (9.*T**2)) )
 return result

#calculate heat capacity and error
#temperature
ePot = energy[:,1]
ePot = ePot[10000::]
T = 2.*eKin.mean()/3.
#calculated by hand from plots
deltaT = 0.0004
#estimate the sample length roughly from Flyvberg&Peterson method
sampleLength = 20000;
heatCap = np.zeros(20)

k =0
for  n in range(0,20):
 sample = ePot[n*sampleLength:(n+1)*sampleLength]
 varSample = np.var(sample)/NUMBER_PARTICLES
 heatCap[n] = heatCapacity(T,varSample)