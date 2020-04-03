#----------------------------------Displacement Driven------------------------------------------#
#--------------------------------------27th March-----------------------------------------------#

#-----------------------------------------------------------------------------------------------#
#                               List of variables used                                          #
#           Ct- Elastic constants                                                               #
#           e - Piezoelectric constants                                                         #
#           k - dielectric constants                                                            #
#-----------------------------------------------------------------------------------------------#

import numpy as np

# A switcher class is used to select the material proprties according to the requirements.
# Three cases are involved 
# Linear Elastic case, Linear Electrical case and Linear Electro-Mechanical case
class Switcher(object):
    def indirect(self,i):
            method_name='LoadCase_'+str(i)
            method=getattr(self,method_name,lambda :'Invalid')
            return method()
    def LoadCase_1(self):
            #---------------------------------------------------------------------#
            #                    Mechanical Case                                  #
            #---------------------------------------------------------------------#

        #----------PZT-PIC151 Material Properties------------#

        #-----------------------Elastic Constants------------------------------#
        Ct = np.array([[110000,64000,0],[64000,100000,0],[0,0,20000]]) #Units MPa

        #---------------Piezoelectric constant matrix-----------------#
        e = np.array([[0,0,0],[0,0,0]])                        #Units C/mm^2

        #---------------Dielectric material constants-----------------#
        k = np.array([[7.54e-12,0],[0,9.82e-12]])                       #Units F/mm
        return Ct,e,k


    def LoadCase_2(self):
            #---------------------------------------------------------------------#
            #                    Electrical Case                                  #
            #---------------------------------------------------------------------#

        #----------PZT-PIC151 Material Properties------------#

        #-----------------------Elastic Constants------------------------------#
        Ct = np.array([[110000,64000,0],[64000,100000,0],[0,0,20000]]) #Units MPa

        #---------------Piezoelectric constant matrix-----------------#
        e = np.array([[0,0,0],[0,0,0]])                        #Units C/mm^2

        #---------------Dielectric material constants-----------------#
        k = np.array([[7.54e-12,0],[0,9.82e-12]])                       #Units F/mm
        return Ct,e,k


    def LoadCase_3(self):
            #---------------------------------------------------------------------#
            #                    Electro Mechanical Case                          #
            #---------------------------------------------------------------------#

        #----------PZT-PIC151 Material Properties------------#

        #-----------------------Elastic Constants------------------------------#
        Ct = np.array([[110000.,64000.,0.],[64000.,100000.,0.],[0.,0.,20000.]]) #Units MPa

        #---------------Piezoelectric constant matrix-----------------#
        e = np.array([[0.,0.,12.0e-6],[-9.6e-6,15.1e-6,0.]])                        #Units C/mm^2

        #---------------Dielectric material constants-----------------#
        k = np.array([[7.54e-12,0],[0,9.82e-12]])                         #Units F/mm
        return Ct,e,k


# Can select the below parameter according to the test case requirement
#   1 - Linear Elastic Loading
#   2 - Linear Electrical Loading
#   3 - Linear Electro-Mechanical Loading
LoadCaseNumber = 3

s=Switcher()
Ct,e,k  = s.indirect(LoadCaseNumber)