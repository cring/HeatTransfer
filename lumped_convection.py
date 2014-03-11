from numpy import *
from matplotlib import pyplot as plt
from scipy import signal

#Geometry
rhoCopper = 8945 # kg/m^3
cCopper = 383.1 #J/kg*K
lCopper = 0.005 #m

areaCopper = 3*lCopper
volumeCopper = lCopper**2

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

print hBlock


def lumped(tempInitial,tempStep,h,n,t):
	temp = []
	time= []

	for i in range(t):
		theta = (tempInitial-tempStep)*e**(-h*n*i)+tempStep
		temp.append([theta])
		time.append([i])
	return temp,time

lumpedBlock = lumped(tempInitial, tempStep, hBlock, nBlock, 100)


plt.xlabel('Time')
plt.ylabel('Temperature')
plt.title('Temperature vs. Time of Lumped Block')
plt.grid(True)
plt.plot(lumpedBlock[1],lumpedBlock[0])
plt.plot([0,100] ,[45,45])
plt.plot([30,30] ,[0,65])
plt.axis([0, 100, 0, 65])
plt.show()