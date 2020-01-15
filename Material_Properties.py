#-----------Before Converting to Displacement Driven-------------------#
#----------- Following code is working for Coupling -------------------#
import numpy as np

yield_stress=210 #MPa
youngs_modulus=210000 #N/mm^2 #MPa
poissons_ratio=0.30

#****************2D Piezoelectric properties******************#

#---------------Piezoelectric constant matrix-----------------#
# e = np.array([[0,0,0],[-6.5e-6,23.3e-6,0]])  #Units C/mm^2
# # Values Taken from Finite Element Modellingand Simulations paper

e = np.array([[0,0,0],[-5.207e-6,15.08e-6,0]])  #Units C/mm^2
# Values taken from Userelement paper

#---------------Dielectric material constants-----------------#
k = np.array([[6.752e-12,0],[0,5.872e-12]])
#k = np.array([[6.752e-12,0],[0,6.752e-12]])
# Values taken from Userelement paper