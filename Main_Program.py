#----------------------------------Displacement Driven------------------------------------------#
#----------------------------One Program for any degree curve-----------------------------------#
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

import sys 

stdoutOrigin=sys.stdout 
sys.stdout = open("log1.txt", "w")

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
        print('Newton Raphson iteration number:',Newton,'\n')
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
            print('Element',j,'nodes',u_e_indices)

            u_e = u_g[u_e_indices]
            #minus 1 since ControlPointAssembly control points numbering starts from 1
            # and python indexing from zero.

            #$$$$$$$$   Defined this in Element routine have to change to main program    $$$$$$$#
            print('Knotconnectivity shape:',np.shape(knotConnectivity))
            print('Knotconnectivity Matrix',knotConnectivity)

            elU = Span_U[knotConnectivity[j,0]]
            elV = Span_V[knotConnectivity[j,1]]
            print('Span range of element:',j,'in xi  direction is:',elU)
            print('Span range of element:',j,'in eta direction is:',elV)

            print('Input Displacement matrix to element Routine:\n',u_e)

            #--------------------Calling Element Routine------------------------------#
            K_e,F_e_int,F_e_ext,sigma,Electrical_Displacement,epsilon,electric_field = elementRoutine(u_e,elU,elV,tau)
            print('Stiffness matrix of element',j,'is:\n',K_e)
            print('Force inetrnal matrix of element',j,'is:\n',F_e_int)
            print('Force external matrix of element',j,'is:\n',F_e_ext)
            print('Stress matrix of element',j,'is:\n',sigma)
            print('Strain matrix of element',j,'is:\n',epsilon)

            # Indices for Force matrix
            # 3 times because of 3 degrees of freedom
            F_e_indices = u_e_indices #Already sorted
            #print('K_e\n',K_e)

            #------------For loops to connect local stiffness matrix to global stiffness matrix------------#
            for val1,index1 in enumerate(u_e_indices):
                print('Connecticity outer loop:',val1,index1)
                for val2,index2 in enumerate(u_e_indices):
                    print('Connecticity inner loop:',val2,index2)
                    Kt_g[index1,index2] = Kt_g[index1,index2] +K_e[val1,val2]
            print('Global Stiffness matrix \n after element no.',j,'is:\n',Kt_g)

########            Kt_g = Kt_g+K_e
            G_global[F_e_indices] = G_global[F_e_indices]+(F_e_int-F_e_ext)
            print('Global Force matrix after element no.',j,'is:\n',G_global)
########            F_g_int = F_g_int+F_e_int
            F_g_int[F_e_indices] = F_g_int[F_e_indices]+F_e_int
            print('Global internal Force matrix after element no.',j,'is:\n',F_g_int)
            print('Fe_Indices:',F_e_indices)

            #$$$$ Previously used Gauss points to plot strain  $$$$$$# 
            #gauss_loc[j] = r_gp
        
        #-----------------Reduced system of equations---------------------#
        K_rg=Kt_g
        print('Global stiffness matrix is:\n',K_rg)
        print('Determinant of Global Stiffness matrix:',np.linalg.det(K_rg))

        #-------------Data from BOunday_Load_Conditions File--------------#
        #-------BCS-----------#
        #BCS:: An Array which is used to delete Global Stiffness matrix rows and coloumns 
        # for solving the equations
        print('Boundary conditions applied are:',BCS)
        K_rg=np.delete(K_rg,BCS,axis=0) #axis=0 is row
        K_rg=np.delete(K_rg,BCS,axis=1) 
        print('Reduced Global stiffness matrix is:\n',K_rg)
        #print('K_rg',K_rg)
        #print('G_global',G_global)
        #print(BCS)
        #print('K_rg',np.linalg.det(K_rg))
        #print('K_rg',K_rg)
        reduced_G_global=G_global
        reduced_G_global=np.delete(reduced_G_global,BCS,axis=0)
        print('Reduced reduced_G_global matrix is:\n',reduced_G_global)
        dU_g=np.matmul(np.linalg.inv(K_rg),-reduced_G_global)
        print('dU_g',dU_g)

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
print('No.of Newton raphson iterations done:',Newton)
print('Reaction force array:\n',F_g_int)
print('Displacements array\n',U_g_0)


#--------------------------Plotting Displacements Calculations----------------------------#
print('control points for the geometry:',P)
print('control points for the geometry P_W:',P_W)
print('P_W[:,:,0:2]',P_W[:,:,0:2])
print('control points for the geometry P_W shape:',np.shape(P_W))
#Seperating out only mechanical displacements for plotting
Mechanical_Displacements=U_g_0[0:ncp*nudof]
#Reshape the mechanical displacements to match with shape of control points
Mechanical_Displacements = Mechanical_Displacements.reshape((ncpeta,ncpxi,2))
print('Mechanical_Displacements',Mechanical_Displacements)
#Adding displacements to control points P_new = P + U  
#P_W is along with corresponding weights
P_W_new_temp = P_W[:,:,0:2] + Mechanical_Displacements
print('P_W_new_temp',P_W_new_temp)
print('P_W_new_temp',P_W_new_temp[:,:,0])
print('P_W_new_temp',P_W_new_temp[:,:,1])
P_W_new = np.zeros_like(P_W)
#P_W_new = np.zeros((ncpeta,ncpxi,4))
P_W_new[:,:,0:2] = P_W_new_temp
P_W_new[:,:,3] = P_W[:,:,3]
print('P_W_new',P_W_new)
P_W_new_T = P_W_new.transpose((1,0,2))
#------------------------------------------------------------------------------------------#


#------------------------Plotting--------------------------#
x_positions=np.linspace(U[0],(U[-1]),ncpxi)
#print('x_positions',x_positions)
y_positions=np.linspace(V[0],(V[-1]),ncpeta)
surface=np.zeros((np.size(x_positions),np.size(y_positions),3))
for i,u in enumerate(x_positions):
    for j,v in enumerate(y_positions):
        # Here send input as displaced positions of control points matrix i.e P_new = P + U
        S_r = NURBS_Surface_Point(n,p,U,m,q,V,P_W_new_T,u,v)
        surface[i,j]=S_r
#Surface generated will be a transpose of what is desired since 
#NURBS_Surface_Point algorithm is modelled in such a way so a transpose is required.
surface = surface.transpose((1,0,2)) 
#print(np.shape(surface))
#print('Surface',surface)
X=surface[:,:,0]
#print('Y',Y)
Y=surface[:,:,1]
print('X',X)
#print('X',X)
#print('Y',Y)
#------------------U1 Plotting----------------------#
fig,ax=plt.subplots()
level = np.linspace(np.amin(X),np.amax(X),13)
plt.contourf(X, Y, X,levels=level,cmap='rainbow')
plt.colorbar()
ax.set(xlabel = 'x [$mm$]', ylabel = 'y [$mm$]', title= 'U1')
#---------------------------------------------------#
#------------------U2 Plotting----------------------#
fig,ax=plt.subplots()
level = np.linspace(np.amin(Y),np.amax(Y),13)
plt.contourf(X, Y, Y,levels=level,cmap='rainbow')
plt.colorbar()
ax.set(xlabel = 'x [$mm$]', ylabel = 'y [$mm$]', title= 'U2')
#---------------------------------------------------#

#------------------RF Plotting Calculations----------------------#
Mechanical_RF = F_g_int[0:ncp*nudof]
Mechanical_RF = Mechanical_RF.reshape((ncpeta,ncpxi,2))
print('Mechanical_RF',Mechanical_RF)
RF1=Mechanical_RF[:,:,0]
print('RF1',RF1)
RF2=Mechanical_RF[:,:,1]
#----------------------------------------------------------------#
#------------------RF1 Plotting----------------------#
fig,ax=plt.subplots()
level = np.linspace(np.amin(RF1),np.amax(RF1),13)
plt.contourf(X, Y, RF1,levels=level,cmap='rainbow')
plt.colorbar()
ax.set(xlabel = 'x [$mm$]', ylabel = 'y [$mm$]', title= 'RF1')
#---------------------------------------------------#

#------------------RF2 Plotting----------------------#
fig,ax=plt.subplots()
level = np.linspace(np.amin(RF2),np.amax(RF2),13)
plt.contourf(X, Y, RF2,levels=level,cmap='rainbow')
plt.colorbar()
ax.set(xlabel = 'x [$mm$]', ylabel = 'y [$mm$]', title= 'RF2')
#---------------------------------------------------#

#-----------------Plotting Electrical Potential Calculations-----------------#
EPOT=U_g_0[nudof*ncp:(nudof*ncp+nedof*ncp)]
EPOT=EPOT.reshape((ncpeta,ncpxi))
print('EPOT',EPOT)
#----------------------------------------------------------------------------#

#-----------------Plotting Electrical Potential -----------------------------#
fig,ax=plt.subplots()
level = np.linspace(np.amin(EPOT),np.amax(EPOT),13)
plt.contourf(X, Y, EPOT,levels=level,cmap='rainbow')
plt.colorbar()
ax.set(xlabel = 'x [$mm$]', ylabel = 'y [$mm$]', title= 'EPOT')
#----------------------------------------------------------------------------#



#--------------Plotting Electrical Reaction Forces Calculations---------------#
Electrical_RF = F_g_int[nudof*ncp:(nudof*ncp+nedof*ncp)]
Electrical_RF = Electrical_RF.reshape((ncpeta,ncpxi))
print('Electrical_RF',Electrical_RF)
#-----------------------------------------------------------------------------#

#-----------------Plotting Electrical Reaction Forces-------------------------#
fig,ax=plt.subplots()
level = np.linspace(np.amin(Electrical_RF),np.amax(Electrical_RF),13)
plt.contourf(X, Y, Electrical_RF,levels=level,cmap='rainbow')
plt.colorbar()
ax.set(xlabel = 'x [$mm$]', ylabel = 'y [$mm$]', title= 'Electrical_RF')
#-----------------------------------------------------------------------------#


sys.stdout.close()
sys.stdout=stdoutOrigin


plt.show()