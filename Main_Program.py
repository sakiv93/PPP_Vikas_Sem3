import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
np.set_printoptions(threshold=np.inf)
from Geometry import *
from Material_Routine import *
from Element_Routine import *
from preprocessing_functions import *

#-----------Import these functions later----------# 
#from mesh_generation import *
#from material_parameters import *
#from assignment_matrix import ASSIGNMENT_MATRIX
#-------------------------------------------------#

#initialization of time like parameter and displacement vector and state variables
tau_start = 0
tau_end = 1
delta_t = 0.1
allTimes = np.arange(tau_start+delta_t, tau_end+delta_t, delta_t)
allTimes[-1] = tau_end
nSteps = len(allTimes)


#Initializing global displcement vector and stress vector 
#necp=Total number of nodes ndof = number of Dof per node 
U_g_0=np.zeros(((nudof+nedof)*necp,1))
global_sigma = np.zeros([nSteps,nel,3,1])


#loop for iterating load from 0% to 100%
#for i in range(nSteps):
for i in range(1):
    #tau = allTimes[i] 
    tau = 1
    u_g = U_g_0 
    #u_g[0,0] = 1/3*E_v*tau*rnodes[0]
    current_sigma = np.zeros([nel,3,1])
    
#------------------DO Newton_Raphson_method----------------------------#
    k=1
    while 1:
        #gauss_loc = np.zeros(nel)
        Kt_g=np.zeros([(nudof+nedof)*necp,(nudof+nedof)*necp])
        G_global=np.zeros(((nudof+nedof)*necp,1))

        # First 8 rows are mechanical forces next 4 rows are electrical forces
        F_g_int=np.zeros(((nudof+nedof)*necp,1))
        F_g_ext=np.zeros(((nudof+nedof)*necp,1))

        
        
#---------------------------Connectivity loop------------------------------#
#
#                           Have to add code
#--------------------------------------------------------------------------#

        for j in range(nel):

            #$$$$$$$$   Defined this in Element routine have to change to main program    $$$$$$$#

            elU = np.array([0,1]) #Manually defined have to generate using connectivity functions
            elV = np.array([0,1])
            #  necp  = number of nodes per element
            u_e = u_g
            #** Call element Routine
            #K_e,F_e_int,F_e_ext, sigma, r_gp = elementRoutine(u_e,tau,rnodes[j:j+2])
            #K_e,F_e_int,F_e_ext, sigma = elementRoutine(u_e,tau)
            K_e,F_e_int, sigma = elementRoutine(u_e,tau)
            #$$$$ Connectivity loop to connect K_e ,F_g_int, F_g_ext to global $$$$#
            Kt_g = Kt_g+K_e
            F_g_int = F_e_int
            #F_g_ext = F_e_ext
            #G_global = G_global+(F_g_int-F_g_ext)
            # current_sigma[j] = sigma 
            # print('Sigma',current_sigma)
            # print('G_global',G_global)
            # print('K_Matrix',Kt_g)

            #$$$$ Previously used Gauss points to plot strain  $$$$$$# 
            #gauss_loc[j] = r_gp
        #Reduced system of equations
        K_rg=Kt_g

        #$$$$  Import BCS File here  $$$$#
        #-------BCS-----------#
        # F_g_ext[5][0] = 100.
        # F_g_ext[7][0] = 100.
        # G_global = G_global+(F_g_int-F_g_ext)
        # u_g[0][0] = 0.
        # u_g[1][0] = 0.
        # u_g[2][0] = 0.
        # u_g[3][0] = 0.
        F_g_ext[2][0] = 100.
        F_g_ext[6][0] = 100.
        G_global = G_global+(F_g_int-F_g_ext)
        u_g[0][0] = 0.
        #u_g[1][0] = 0.
        u_g[4][0] = 0.
        u_g[5][0] = 0.
        # Electrical bcs
        u_g[8][0] = 0.
        u_g[10][0] = 0.
        K_rg=np.delete(K_rg,[0,4,5,8,10],axis=0) #axis=0 is row
        K_rg=np.delete(K_rg,[0,4,5,8,10],axis=1)
        # K_rg=np.delete(K_rg,[0,1,2,3],axis=0) #axis=0 is row
        # K_rg=np.delete(K_rg,[0,1,2,3],axis=1)
        reduced_G_global=G_global
        reduced_G_global=np.delete(reduced_G_global,[0,4,5,8,10],axis=0)
        # reduced_G_global=np.delete(reduced_G_global,[0,1,2,3],axis=0)
        dU_g=np.matmul(np.linalg.inv(K_rg),-reduced_G_global)
        #dU_g_insert=np.insert(dU_g,[0,1,2,3],0,axis=0)
        # dU_g_insert=np.insert(dU_g,[0],0,axis=0)
        # dU_g_insert=np.insert(dU_g_insert,[0],0,axis=0)
        # dU_g_insert=np.insert(dU_g_insert,[0],0,axis=0)
        # dU_g_insert=np.insert(dU_g_insert,[0],0,axis=0)
        dU_g_insert=np.insert(dU_g,[0],0,axis=0)
        #dU_g_insert=np.insert(dU_g_insert,[1],0,axis=0)
        dU_g_insert=np.insert(dU_g_insert,[4],0,axis=0)
        dU_g_insert=np.insert(dU_g_insert,[5],0,axis=0)
        dU_g_insert=np.insert(dU_g_insert,[8],0,axis=0)
        dU_g_insert=np.insert(dU_g_insert,[10],0,axis=0)
        u_g=u_g+dU_g_insert
        # print('F_External',F_g_ext)
        # print('F_Internal',F_g_int)
        # print('U Global',u_g)
        # print('Reduced K matrix',K_rg)

        if (np.linalg.norm(reduced_G_global,np.inf)<0.005*np.linalg.norm(F_g_int,np.inf) or np.linalg.norm(dU_g,np.inf)<0.005*np.linalg.norm(u_g[1:],np.inf)) or k > 5 :     
            break
        else:
            k=k+1
    if k>5:
        print("Convergence criterion is not met, at load (peercentage):",tau*100)
        break
    else:
        U_g_0 = u_g
        global_sigma[i] = current_sigma
print(U_g_0)
#print(k)
