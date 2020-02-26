#----------------------------------Displacement Driven------------------------------------------#
#------------------------------------------26 Feb-----------------------------------------------#
import numpy as np

yield_stress=210 #MPa
youngs_modulus=210000 #N/mm^2 #MPa
poissons_ratio=0.30

#****************2D Piezoelectric properties******************#

#---------------Piezoelectric constant matrix-----------------#
# e = np.array([[0,0,0],[-6.5e-6,23.3e-6,0]])  #Units C/mm^2
# # Values Taken from Finite Element Modellingand Simulations paper

#e = np.array([[0,0,0],[0,0,0]])
e = np.array([[0,0,1.271e-5],[-5.207e-6,15.08e-6,0]])  #Units C/mm^2
# Values taken from Userelement paper

#---------------Dielectric material constants-----------------#
k = np.array([[6.752e-12,0],[0,5.872e-12]])
#k = np.array([[0,0],[0,0]])
#k = np.array([[6.752e-12,0],[0,6.752e-12]])
# Values taken from Userelement paper








# #----------PZT-PIC151 Material Properties------------#
# #-----------------------Elastic Constants------------------------------#
# yield_stress=210 #MPa
# youngs_modulus=210000 #N/mm^2 #MPa
# poissons_ratio=0.30
# #Ct = np.array([[110000,64000,0],[64000,100000,0],[0,0,20000]]) #Units MPa

# #---------------Piezoelectric constant matrix-----------------#
# e = np.array([[0,0,12.0],[-9.6,15.1,0]])  
# #e = np.array([[0,0,0],[0,0,0]])                   #Units C/m^2

# #---------------Dielectric material constants-----------------#
# k = np.array([[0.00754,0],[0,0.00982]]) 
# #k = np.array([[0,0],[0,0]])                     #Units microF/m