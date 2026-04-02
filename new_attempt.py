import numpy as np
import math ## needed this for code to work on my version
import matplotlib.pyplot as plt
import re

#****************************** Disk Parameters****************************************
 
#Define the power law, scaling for the disk emission, and number of radial steps from .07 AU to 200. AU, with step size .05 AU
tpow=-0.5	#defines the power law for the disk temperature profile
Dscale=0.35   #scales the strength of the disk emission relative to the star
Rstepsize=0.05   #the width of each annulus in AU
Rsub=0.07   #the sublimation radius
Rend=300  #ending the disk at 300AU
RstepNumber= int(round( (Rend-Rsub)/Rstepsize  ) ) #Gives integer number of steps from 0.7-200 AU
 
#Define numpy arrays for starting radius, area, and temperature of each annulus, then fill the arrays
R=[] # range(RstepNumber)  # will be inner radius of annulus
Area=[] #=range(RstepNumber-1) # area of annulus , since total lum  ~ I x area
Temp=[] #=range(RstepNumber-1) #Average temperature of our annuli

Rdisk = 4
for j in range(0,RstepNumber):
  R.append(Rdisk)
  Rdisk= Rdisk + Rstepsize
 
for j in range(0,RstepNumber-1):
  Area.append((R[j+1]**2- R[j]**2))
  Rmid= (R[j+1] + R[j])/2.0
  Temp.append(1800*(Rmid/Rsub)**(tpow))

  #***** FLUX FROM DISK *************
 
wstart=10e-9  # starting wavelength in m
wend=1000e-6   #end wavelength in m
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

# DATA POINTS FOR GM AUR, wavelength in microns, lambda Flambda in ergs/s/cm^2
 
GMAur_L=[0.36890015,0.45258972,5.86e-01,0.666506,1.3083931,1.7492352,2.2092927,3.6587987,4.4578514,5.849414,12.335447,27.036976,62.506542,107.64069,631.3774,1079.9856]
GMAur_F=[7.08e-11,1.75e-10,3.36e-10,5.01e-10,6.91e-10,6.08e-10,3.32e-10,1.24e-10,8.89e-11,4.87e-11,6.36e-11,1.10e-10,1.47e-10,9.55e-11,6.20e-12,1.01e-12]
print(GMAur_F)
print(GMAur_L)

GMAur_Ffixed = np.array(GMAur_F, dtype=np.float64)

## Flux (F) at given wavelengths (L)
##Make the flux match at relevant wavelengths
##you have to match it to the <1 microns and say that that is ONLY the sun's impact.

####### My attempt at repeating the process for GMAur #######
##First, we must normalize!!

##C = Ftotal/GMAur_F
##I'm going to choose the first GMAur_L point for my ~.35 micron point, and the first GMAur_F point for my relative flux at that point.
C = Ftotal[3]/GMAur_F[0]

#print(wavePlot[35]) ## this is the array i found initially. it's 10x too big
print(wavePlot[3]) ## so the relative flux is Ftotal[3]


## i played around with the range here, frmo printing the whole array, to limiting it lower and higher.
#for i in range(0, 15):
#  print(wavePlot[i])


GMAur_F= C*np.array(GMAur_F) #Scale the GM Aur fluxes to match model SED fluxes at ~0.35 micron

##Back to lab notes

# Change wavelengths to microns, and define new arrays that equal "lambdaxF_lambda"
wavePlot=np.array(wavePlot)*10**6
LFdisk=np.array(Fdisk)*np.array(wavePlot)
LFstar=np.array(Fstar)*np.array(wavePlot)
LFtotal=np.array(Ftotal)*np.array(wavePlot)
 
plt.plot(wavePlot,np.log10(LFdisk), label="Flared Disk")
plt.plot(wavePlot,np.log10(LFstar), label="Star")
plt.plot(wavePlot,np.log10(LFtotal), label="Disk Plus Star")
plt.plot(GMAur_L,np.log10(GMAur_F), label="GMAuriga Flux")
plt.xlim([.1,2000])
plt.xscale('log')
plt.ylim([-10,0])
plt.legend(loc='best',prop={'size': 10})
plt.xlabel('Wavelength (λ)')
plt.ylabel('Flux Normalized with Lambda')
plt.show()
