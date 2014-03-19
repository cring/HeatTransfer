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

lumpedBlock = lumped(tempInitial, tempStep, hBlock, nBlock, 100)

def coupledLumped(tempInitial,tempStep,h,t):
	C1 = h*areaCopper/(rhoCopper*cCopper*volumeCopper)
	C2 = kAcrylic*lCopper/(sqrt(alphaAcrylic)*rhoCopper*cCopper*volumeCopper)

	D1 = C2/2-cmath.sqrt(C2**2-4*C1)/2
	D2 = C2/2+cmath.sqrt(C2**2-4*C1)/2

	print D1
	print s.erfc(D1)
	print e**(D1**2)
	print C1/((D2-D1)*(D2*D1))
	print D1*e**(D2**2*1)*s.erfc(D2*sqrt(1))
	print (D2*e**(D1**2*1)*s.erfc(D1*sqrt(1)))
	print C1/((D2-D1)*(D2*D1))*((D1*e**(D2**2*1)*s.erfc(D2*sqrt(1)))-(D2*e**(D1**2*1)*s.erfc(D1*sqrt(1)))+D2-D1)

	temp = []
	time= []

	for i in range(t):
		temp_t = (tempStep-tempInitial)*C1/((D2-D1)*(D2*D1))*((D1*e**(D2**2*i)*s.erfc(D2*sqrt(i)))-(D2*e**(D1**2*i)*s.erfc(D1*sqrt(i)))+D2-D1)+tempInitial
		temp.append([temp_t.real])
		time.append([i])
	return temp,time

coupledlumpedBlock = coupledLumped(tempInitial, tempStep, 153.8, 100)

print hBlock
print coupledlumpedBlock[0][30]


plt.xlabel('Time')
plt.ylabel('Temperature')
plt.title('Temperature vs. Time of Lumped Block')
plt.grid(True)
plt.plot(lumpedBlock[1],lumpedBlock[0])
plt.plot(coupledlumpedBlock[1],coupledlumpedBlock[0])
plt.plot([0,100] ,[45,45])
plt.plot([30,30] ,[0,65])
plt.axis([0, 100, 0, 65])
plt.show()