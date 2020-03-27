#----------------------------------Displacement Driven------------------------------------------#
#--------------------------------------27th March-----------------------------------------------#
import numpy as np

class Switcher(object):
    def indirect(self,i):
            method_name='LoadCase_'+str(i)
            method=getattr(self,method_name,lambda :'Invalid')
            return method()
    def LoadCase_1(self):
            #---------------------------------------------------------------------#
            #                    Pure Mechanical Case                             #
            #---------------------------------------------------------------------#

        #----------PZT-PIC151 Material Properties------------#

        #-----------------------Elastic Constants------------------------------#
        Ct = np.array([[110000,64000,0],[64000,100000,0],[0,0,20000]]) #Units MPa

        #---------------Piezoelectric constant matrix-----------------#
        e = np.array([[0,0,0],[0,0,0]])                        #Units C/m^2

        #---------------Dielectric material constants-----------------#
        k = np.array([[7.54e-12,0],[0,9.82e-12]])                       #Units microF/m
        return Ct,e,k


    def LoadCase_2(self):
            #---------------------------------------------------------------------#
            #                    Pure Electrical Case                             #
            #---------------------------------------------------------------------#

        #----------PZT-PIC151 Material Properties------------#

        #-----------------------Elastic Constants------------------------------#
        Ct = np.array([[110000,64000,0],[64000,100000,0],[0,0,20000]]) #Units MPa

        #---------------Piezoelectric constant matrix-----------------#
        e = np.array([[0,0,0],[0,0,0]])                        #Units C/m^2

        #---------------Dielectric material constants-----------------#
        k = np.array([[7.54e-12,0],[0,9.82e-12]])                       #Units microF/m
        return Ct,e,k


    def LoadCase_3(self):
            #---------------------------------------------------------------------#
            #                    Electro Mechanical Case                          #
            #---------------------------------------------------------------------#

        #----------PZT-PIC151 Material Properties------------#

        #-----------------------Elastic Constants------------------------------#
        Ct = np.array([[110000.,64000.,0.],[64000.,100000.,0.],[0.,0.,20000.]]) #Units MPa

        #---------------Piezoelectric constant matrix-----------------#
        e = np.array([[0.,0.,12.0e-6],[-9.6e-6,15.1e-6,0.]])                        #Units C/m^2

        #---------------Dielectric material constants-----------------#
        k = np.array([[7.54e-12,0],[0,9.82e-12]])                         #Units microF/m
        return Ct,e,k


LoadCaseNumber = 3


s=Switcher()
Ct,e,k  = s.indirect(LoadCaseNumber)
print(Ct)
print(e)
print(k)



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


# #---------------------------------------------------------------------#
# #                    Electro Mechanical Case                          #
# #---------------------------------------------------------------------#

# yield_stress=210 #MPa
# youngs_modulus=210000 #N/mm^2 #MPa
# poissons_ratio=0.30
# #---------------Piezoelectric constant matrix-----------------#
# # Give here reference for values taken from
# e = np.array([[0,0,1.271e-5],[-5.207e-6,15.08e-6,0]])  #Units C/mm^2
# #---------------Dielectric material constants-----------------#
# # Give here reference for values taken from
# k = np.array([[6.752e-12,0],[0,5.872e-12]])             #Units Write units
# return yield_stress,youngs_modulus,poissons_ratio,e,k