#-------------------------------Displacement Driven----------------------------------------#
#----------- Following code is working for Displacement Driven Coupling -------------------#
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
np.set_printoptions(threshold=np.inf)
from Geometry import *
from Material_Routine import *
from Element_Routine import *
from preprocessing_functions import *
from Boundary_Load_Conditions import *

#initialization of time like parameter and displacement vector and state variables
tau_start = 0
tau_end = 1
delta_t = 0.1
allTimes = np.arange(tau_start+delta_t, tau_end+delta_t, delta_t)
allTimes[-1] = tau_end
nSteps = len(allTimes)

#Initializing global Stress Tensor for values of stress at gauss points in each step
#***** Have to Store values at each gauss point ******#
#global_sigma = np.zeros([nSteps,nel,3,1])


#loop for iterating load from 0% to 100%
# for i in range(nSteps):
for i in range(1):
    #tau = allTimes[i] 
    tau = 1
    u_g = U_g_0
    #***** Have to Store values at each gauss point ******#
    #current_sigma = np.zeros([nel,3,1])
    
#------------------DO Newton_Raphson_method----------------------------#
    Newton=1
    while 1:
        #$$$$ Previously used Gauss points to plot strain  $$$$$$# 
        #gauss_loc = np.zeros(nElem)
        Kt_g=np.zeros([(nudof+nedof)*necp,(nudof+nedof)*necp])
        G_global=np.zeros(((nudof+nedof)*necp,1))

        # First 8 rows are mechanical forces next 4 rows are electrical forces
        F_g_int=np.zeros(((nudof+nedof)*necp,1))

#---------------------------Connectivity loop------------------------------#
#
#                           Have to add code
#--------------------------------------------------------------------------#

        for j in range(nel):

            #$$$$$$$$   Defined this in Element routine have to change to main program    $$$$$$$#

            elU = np.array([0,1]) #Manually defined have to generate using connectivity functions
            elV = np.array([0,1])
            u_e = u_g
            #print('Input Displacement matrix to element Routine:',u_e)

            #--------------------Calling Element Routine------------------------------#
            K_e,F_e_int,F_e_ext,sigma,Electrical_Displacement,epsilon = elementRoutine(u_e,tau)
            #print(K_e)
            #print('u_e output from element routine:',u_e)
            #$$$$ Connectivity loop to connect K_e ,F_g_int, F_g_ext to global $$$$#
            Kt_g = Kt_g+K_e
            G_global = G_global+(F_e_int-F_e_ext)
            F_g_int = F_g_int+F_e_int

            #$$$$ Previously used Gauss points to plot strain  $$$$$$# 
            #gauss_loc[j] = r_gp
        
        #-----------------Reduced system of equations---------------------#
        K_rg=Kt_g

        #-------------Data from BOunday_Load_Conditions File--------------#
        #-------BCS-----------#
        #BCS:: An Array which is used to delete Global Stiffness matrix rows and coloumns 
        # for solving the equations
        K_rg=np.delete(K_rg,BCS,axis=0) #axis=0 is row
        K_rg=np.delete(K_rg,BCS,axis=1) 
        #print('K_rg',K_rg)
        #print('G_global',G_global)
        reduced_G_global=G_global
        reduced_G_global=np.delete(reduced_G_global,BCS,axis=0)
        dU_g=np.matmul(np.linalg.inv(K_rg),-reduced_G_global)
        #print('dU_g',dU_g)

        #-------For Newton Raphson Scheme covergence Criterion--------------#
        dU_g_convergence=dU_g
        u_g_convergence=np.delete(u_g,BCS,axis=0)
        dU_g_insert=dU_g
        #----dU_g_insert matrix is inserted with Boundary and Load conditions nodal values as 0------#
        for val in BCS:
            dU_g_insert=np.insert(dU_g_insert,val,0,axis=0)
            #print(dU_g_insert)
        #print('dU_g_insert',dU_g_insert)
        u_g=u_g+dU_g_insert

        if (np.linalg.norm(reduced_G_global,np.inf)<0.005*np.linalg.norm(F_g_int,np.inf) or np.linalg.norm(dU_g_convergence,np.inf)<0.005*np.linalg.norm(u_g_convergence,np.inf) or Newton > 5) :     
            break
        else:
            Newton=Newton+1
    if Newton>5:
        print("Convergence criterion is not met, at load (peercentage):",tau*100)
        break
    else:
        U_g_0 = u_g
        #print('U For next Iteration:',U_g_0)
        #global_sigma[i] = current_sigma
print('F_internal',F_g_int)
print('Displacements',U_g_0)
print(P)
print('Sigma',sigma)
print('Strain',epsilon)
