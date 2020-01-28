#----------------------------------Displacement Driven------------------------------------------#
#----------- Following code is working for Coupling with connectivity matrix in progress-------------------#
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
        file=open('Random.txt','a+')
        file.write('Newton Raphson iteration number:%d \n' % (Newton))
        file.close()
        #$$$$ Previously used Gauss points to plot strain  $$$$$$# 
        #gauss_loc = np.zeros(nElem)
########        Kt_g=np.zeros([(nudof+nedof)*necp,(nudof+nedof)*necp])
        Kt_g=np.zeros([(nudof+nedof)*ncp,(nudof+nedof)*ncp])
########        G_global=np.zeros(((nudof+nedof)*necp,1))
        G_global=np.zeros(((nudof+nedof)*ncp,1))

        # First 8 rows are mechanical forces next 4 rows are electrical forces
########        F_g_int=np.zeros(((nudof+nedof)*necp,1))
        F_g_int=np.zeros(((nudof+nedof)*ncp,1))

#---------------------------Connectivity loop------------------------------#
#
#                           Have to add code
#--------------------------------------------------------------------------#

        for j in range(nel):
            file=open('Random.txt','a+')
            file.write('Element number:%d \n' % (j+1))
            file.close()

            #----- If p=1 or q=1 -------#
            u_e_indices=ControlPointAssembly(ncpxi,p,ncpeta,q,j)
            # #----- If p>1 or q>1 -------#
            # u_e_indices=ControlPointAssembly(ncpxi,p_ord,ncpeta,q_ord,nel)
            #print('u_e_indices',u_e_indices)


            # Generated values of 'u_e_indices' are of float. 
            # Converting to int for using for using it as indices of array
            u_e_indices = u_e_indices.astype(int)
            #u_e_indices = np.concatenate((2*u_e_indices+0,2*u_e_indices+1,(ncp*2)*u_e_indices))  
            u_e_indices = np.concatenate((nudof*u_e_indices+0,nudof*u_e_indices+1,(ncp*nudof)+u_e_indices))
            u_e_indices=np.sort(u_e_indices)
            #print('u_e_indices',u_e_indices)

            with open('Random.txt', 'a+') as f:
                f.write('Element %d indices:' % (j))
                for item in u_e_indices:
                    f.write("%s," % item)
                f.write('\n')

            u_e = u_g[u_e_indices]
            #minus 1 since ControlPointAssembly control points numbering starts from 1
            # and python indexing from zero.

            #$$$$$$$$   Defined this in Element routine have to change to main program    $$$$$$$#
            #print('Knotconnectivity shape',np.shape(knotConnectivity))
            #print(knotConnectivity)
            with open('Random.txt', 'a+') as f:
                f.write('Knotconnectivity:')
                for item in knotConnectivity[j-1]:
                    f.write("%s," % item)
                f.write('\n')
            print('Span_U:',Span_U)
            elU = Span_U[knotConnectivity[j,0]]
            elV = Span_V[knotConnectivity[j,1]]
            with open('Random.txt', 'a+') as f:
                f.write('Span_U :')
                for item in Span_U[knotConnectivity[j-1,0]]:
                    f.write("%s," % item)
                f.write('\n')

            with open('Random.txt', 'a+') as f:
                f.write('Span_V :')
                for item in Span_V[knotConnectivity[j-1,1]]:
                    f.write("%s," % item)
                f.write('\n')

########            elU = np.array([0,1]) #Manually defined have to generate using connectivity functions
########            elV = np.array([0,1])
########            u_e = u_g
            #print('Input Displacement matrix to element Routine:',u_e)

            #--------------------Calling Element Routine------------------------------#
            K_e,F_e_int,F_e_ext,sigma,Electrical_Displacement,epsilon = elementRoutine(u_e,elU,elV,tau)
            #print(K_e)
            #print('u_e output from element routine:',u_e)
            #$$$$ Connectivity loop to connect K_e ,F_g_int, F_g_ext to global $$$$#

            # Indices for Force matrix
            # 3 times because of 3 degrees of freedom
            F_e_indices = u_e_indices #Already sorted
            #print('K_e\n',K_e)

            #------------For loops to connect local stiffness matrix to global stiffness matrix------------#
            for val1,index1 in enumerate(u_e_indices):
                for val2,index2 in enumerate(u_e_indices):
                    Kt_g[index1,index2] = Kt_g[index1,index2] +K_e[val1,val2]
            #print('Kt_g\n',Kt_g)

########            Kt_g = Kt_g+K_e
            G_global[F_e_indices] = G_global[F_e_indices]+(F_e_int-F_e_ext)
########            F_g_int = F_g_int+F_e_int
            F_g_int[F_e_indices] = F_g_int[F_e_indices]+F_e_int

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
        #print(BCS)
        #print('K_rg',np.linalg.det(K_rg))
        #print('K_rg',K_rg)
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
print('Displacements',U_g_0[[0,1,4,5,12,13,24,16,17,18,20,24,26]])
