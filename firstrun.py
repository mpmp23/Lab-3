import numpy as np
import math ## needed this for code to work on my version
import matplotlib.pyplot as plt
import re

#****************************** Disk Parameters****************************************
 
#Define the power law, scaling for the disk emission, and number of radial steps from .07 AU to 200. AU, with step size .05 AU
tpow=-0.75	#defines the power law for the disk temperature profile
Dscale=0.25   #scales the strength of the disk emission relative to the star
Rstepsize=0.05   #the width of each annulus in AU
Rsub=0.07    #the sublimation radius
Rend=200  #ending the disk at 200AU
RstepNumber= int(round( (Rend-Rsub)/Rstepsize  ) ) #Gives integer number of steps from 0.7-200 AU
 
#Define numpy arrays for starting radius, area, and temperature of each annulus, then fill the arrays
R=[] # range(RstepNumber)  # will be inner radius of annulus
Area=[] #=range(RstepNumber-1) # area of annulus , since total lum  ~ I x area
Temp=[] #=range(RstepNumber-1) #Average temperature of our annuli

Rdisk = Rsub
for j in range(0,RstepNumber):
  R.append(Rdisk)
  Rdisk= Rdisk + Rstepsize
 
for j in range(0,RstepNumber-1):
  Area.append((R[j+1]**2- R[j]**2))
  Rmid= (R[j+1] + R[j])/2.0
  Temp.append(1800*(Rmid/Rsub)**(tpow))

  #***** FLUX FROM DISK *************
 
wstart=10e-9  # starting wavelength in m
wend=500e-6   #end wavelength in m
wstep =100e-9  #step for wavelength in m
wstepNumber=int(round( (wend-wstart)/wstep) )
 
wavePlot=[] # save wavelengths in array for easy plotting
Fdisk=[]  # will be the total disk luminosity per wavelength, summed up from all the annuli
Fstar=[]  # the star's luminosity per wavelength
 
Sumstar=0
Sumdisk=0

#define some values for the Planck function
 
def PlanckLaw(T, wave):
    c=1.1927*10**-16 #2hc^2
    d=.014394135  #hc/Kb
    if (d/(T*wave))>500:  #Put in this condition to avoid getting overflow of digits in exponential
        return 0
    else:
        B= ((c/wave)**5)*(math.exp(d/(T*wave)) - 1 )**-1
        return B

#Fill the wavelength, disk luminosity, and star luminosity arrays
wave=wstart
for i in range(0, wstepNumber):
    wavePlot.append(wave)
    F_annuli=0
    Fstar.append(PlanckLaw(4000,wavePlot[i]))
    Sumstar=Fstar[i] + Sumstar
    for j in range(0,RstepNumber-1):
            F_annuli=PlanckLaw(Temp[j],wavePlot[i])*Area[j] + F_annuli

    Fdisk.append(F_annuli)
    wave= wave + wstep
    Sumdisk=Fdisk[i]+Sumdisk

 
# NORMALIZE ALL THE FLUXES and scale disk flux strength to DScale
Fstar=np.array(Fstar)/Sumstar
Fdisk=Dscale*np.array(Fdisk)/Sumdisk
Ftotal= Fstar + Fdisk

# Change wavelengths to microns, and define new arrays that equal "lambdaxF_lambda"
wavePlot=np.array(wavePlot)*10**6
LFdisk=np.array(Fdisk)*np.array(wavePlot)
LFstar=np.array(Fstar)*np.array(wavePlot)
LFtotal=np.array(Ftotal)*np.array(wavePlot)
 
plt.plot(wavePlot,np.log10(LFdisk), label="flat disk")
plt.plot(wavePlot,np.log10(LFstar), label="star")
plt.plot(wavePlot,np.log10(LFtotal), label="disk plus star")
plt.xlim([.1,2000])
plt.xscale('log')
plt.ylim([-10,0])
plt.legend(loc='best',prop={'size': 10})
plt.xlabel('wavelength (lambda)')
plt.ylabel('Flux normalized with lambda')
plt.show()