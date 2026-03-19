

# DATA POINTS FOR GM AUR, wavelength in microns, lambda Flambda in ergs/s/cm^2
 
GMAur_L=[0.36890015, 0.45258972  ,5.86E-01  ,  0.666506  ,  1.3083931  ,  1.7492352 ,   2.2092927 ,   3.6587987  , 4.4578514, 5.849414,  12.335447  ,  27.036976  ,  62.506542  ,  107.64069 ,   631.3774 ,   1079.9856]
GMAur_F=[7.08E-11, 1.75E-10, 3.36E-10 , 5.01E-10 , 6.91E-10 , 6.08E-10 , 3.32E-10 ,1.24E-10 , 8.89E-11, 4.87E-11, 6.36E-11 , 1.10E-10  , 1.47E-10  , 9.55E-11  , 6.20E-12  , 1.01E-12 ]
 
GMAur_F= *np.array(GMAur_F) #Scale the GM Aur fluxes to match model SED fluxes at ~0.35 micron
