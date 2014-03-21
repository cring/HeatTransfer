from numpy import *
from matplotlib import pyplot as plt
from scipy import signal
from scipy import special as s
import cmath

#Geometry Copper
rhoCopper = 8945 # kg/m^3
cCopper = 383.1 #J/kg*K
lCopper = 0.005 #m

areaCopper = 3*lCopper
volumeCopper = lCopper**2

#Properties Acrylics -- http://www.builditsolar.com/References/Glazing/physicalpropertiesAcrylic.pdf
kAcrylic = 0.19 #W/mK
cAcrylic = 1470 #J/kg*k
rhoAcrylic = 1.19*1000 #kg/m^3
alphaAcrylic = kAcrylic/(rhoAcrylic*cAcrylic)

#Initial Temperature
tempInitial = 25.0 #C
tempStep = 65.0 #C
tempBlock = 45.0 #C

#Temp Differences
thetaCurrent = tempBlock-tempStep
thetaInitial = tempInitial-tempStep

#Lumped Capacitance Theta(t)/Theta0 = exp(-n*h*t)
nBlock = areaCopper/(rhoCopper*volumeCopper*cCopper)

#Calculate h at time = 3
currentTime = 30
hBlock = -log(thetaCurrent/thetaInitial)/(nBlock*currentTime)

def lumped(tempInitial,tempStep,h,n,t):
	temp = []
	time= []

	for i in range(t):
		temp_t = (tempInitial-tempStep)*e**(-h*n*i)+tempStep
		temp.append([temp_t])
		time.append([i])
	return temp,time

def coupledLumped(tempInitial,tempStep,h,t):
	C1 = h*areaCopper/(rhoCopper*cCopper*volumeCopper)
	C2 = kAcrylic*lCopper/(sqrt(alphaAcrylic)*rhoCopper*cCopper*volumeCopper)

	D1 = C2/2-cmath.sqrt(C2**2-4*C1)/2
	D2 = C2/2+cmath.sqrt(C2**2-4*C1)/2

	temp_t = (tempStep-tempInitial)*C1/((D2-D1)*(D2*D1))*((D1*e**(D2**2*t)*s.erfc(D2*sqrt(t)))-(D2*e**(D1**2*t)*s.erfc(D1*sqrt(t)))+D2-D1)+tempInitial
	
	return temp_t.real
	
def coupledCalculator(tempInitial,tempStep,h,t):
	calcTemp=[]
	calcTime=[]
	for i in range(t):		
		calcTemp.append(coupledLumped(tempInitial, tempStep, h, i))
		calcTime.append([i])
	return calcTemp,calcTime

def timeDelayCalculator(tempInitial,tempStep,h,start,stop,step):
	calcTemp = []
	calcTime = []

	for i in arange(0,stop,step):
		if i < start:
			calcTemp.append(0)
			calcTime.append([i])

		else:
			calcTemp.append(coupledLumped(tempInitial, tempStep, h, i-start))
			calcTime.append([i])
	
	return calcTemp,calcTime

def timeSumCalculator():
	timeStepList = [0,0.2,0.4,0.6,0.8,1]
	tempStepList = [25,43,55,63,65,65]

	summedCalc = zeros(1000).tolist()

	for i in range(len(timeStepList)-1):
		tempInitial = 0
		tempStep = tempStepList[i+1]-tempStepList[i]
		convolution = timeDelayCalculator(tempInitial,tempStep,153.8,timeStepList[i],100,0.1)

		for i in range(len(convolution[0])):
			summedCalc[i] += convolution[0][i]

	for i in range(len(summedCalc)):
			summedCalc[i] += tempStepList[0]
	
	return summedCalc , convolution[1]


lumpedBlock = lumped(tempInitial, tempStep, hBlock, nBlock, 100)

#coupledlumpedBlock2 = coupledCalculator(tempInitial, tempStep, hBlock, 100)
coupledlumpedBlock2 = coupledCalculator(tempInitial, tempStep, 153.8, 100)
coupledlumpedBlock = timeSumCalculator()

#print coupledlumpedBlock
print len(coupledlumpedBlock[0])
print len(coupledlumpedBlock[1])

plt.xlabel('Time')
plt.ylabel('Temperature')
plt.title('Temperature vs. Time of Lumped Block')
plt.grid(True)
plt.plot(lumpedBlock[1],lumpedBlock[0],label = "Lumped - h~131")
plt.plot(coupledlumpedBlock[1],coupledlumpedBlock[0], label = "Coupled Lumped - h~154")
plt.plot(coupledlumpedBlock2[1],coupledlumpedBlock2[0], label = "Uncorrected Coupled Lumped - h~131")
#plt.legend(loc=4)
plt.plot([0,100] ,[45,45])
plt.plot([30,30] ,[0,65])
plt.axis([0, 100, 20, 65])
plt.show()