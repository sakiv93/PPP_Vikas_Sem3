#----------------------------------Displacement Driven------------------------------------------#
#--------------------------------------27th March-----------------------------------------------#

#-----------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------#
#                               List of variables used                                                      #
#           p,q                     - Degree of the curve in xi and eta direction respectively              #
#           xi,eta                  - Parametric coordinates in xi and eta direction respectively           #
#           ncpxi, ncpeta           - No.of control points in xi and eta direction respectively             #
#           Pw                      - Control points matrix                                                 #
#           Pw_new                  - Deformed Control points matrix (P+U)                                  # 
#           XI, ETA                 - Knot Vectors in xi and eta direction respectively                     #
#           ncp                     - Total number of control points                                        #
#           necp                    - Total number of control points in an element                          # 
#           nel                     - Total number of elements                                              #
#           sigma                   - stress                                                                #
#           Electrical_Displacement - Electrical Displacement                                               #
#           epsilon                 - strain values                                                         #
#           electric_field          - electric field values                                                 #
#           electric_field          - electric field values                                                 #
#           U_g_0                   - Stores final value of DOF's                                           #
#           u_e                     - Stores values of element DOF's                                        #
#           kt_g                    - Global  stiffness matrix                                              #
#           kt_rg                   - Reduced Global  stiffness matrix                                      #
#           K_e                     - Element stiffness matrix                                              #
#           F_g_int                 - Global  internal force matrix                                         #
#           F_e_int                 - Element internal force matrix                                         #
#           F_e_ext                 - Element external force matrix                                         #
#           u_e_indices             - Stores the element DOF numbers                                        #
#-----------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt
from Geometry import *
from Material_Routine import *
from Element_Routine import *
from preprocessing_functions import *
from Boundary_Load_Conditions import *
from pathlib import Path
import sys 
import timeit

stdoutOrigin=sys.stdout 
sys.stdout = open("log.txt", "w")

Main_start = timeit.timeit()
for i in range(1):
    tau = 1
    u_g = U_g_0
    
    #------------------DO Newton_Raphson_method----------------------------#
    Newton=1
    while 1:
        #-----------Initializing------------# 
        Kt_g=np.zeros([(nudof+nedof)*ncp,(nudof+nedof)*ncp])    # Global Stiffness matrix 
        G_global=np.zeros(((nudof+nedof)*ncp,1))                # Global G (G = F_int-F_ext)

        #First 2*ncp rows are mechanical forces next 4 rows are electrical forces in each element 
        F_g_int=np.zeros(((nudof+nedof)*ncp,1))                 # Global Force internal

        #--------Looping over elements in the simulation----------# 
        for j in range(nel):
            # Fetching indices of displacement DOF in an element.
            u_e_indices=ControlPointAssembly(ncpxi,p,ncpeta,q,j)

            # Generated values of 'u_e_indices' are of float. 
            # Converting to int for using it as indices of array
            u_e_indices = u_e_indices.astype(int)
            u_e_indices = np.concatenate((nudof*u_e_indices+0,nudof*u_e_indices+1,(ncp*nudof)+u_e_indices))
            u_e_indices=np.sort(u_e_indices)
            u_e = u_g[u_e_indices]

            #------ Extracting values of the knot span for the element in each direction
            elXI =  Span_XI[knotConnectivity[j,0]]
            elETA = Span_ETA[knotConnectivity[j,1]]

            #--------------------Calling Element Routine------------------------------#
            K_e,F_e_int,F_e_ext,sigma,Electrical_Displacement,epsilon,electric_field = elementRoutine(u_e,elXI,elETA,tau)

            F_e_indices = u_e_indices # Already sorted

            #------------For loops to connect local stiffness matrix to global stiffness matrix------------#
            for val1,index1 in enumerate(u_e_indices):
                #print('Connecticity outer loop:',val1,index1)
                for val2,index2 in enumerate(u_e_indices):
                    #print('Connecticity inner loop:',val2,index2)
                    Kt_g[index1,index2] = Kt_g[index1,index2] +K_e[val1,val2]
            G_global[F_e_indices] = G_global[F_e_indices]+(F_e_int-F_e_ext)
            F_g_int[F_e_indices] = F_g_int[F_e_indices]+F_e_int
        ElementLoop_end = timeit.timeit()
        #--------End of Looping over elements in the simulation----------#
        
        #-----------------Reduced system of equations---------------------#
        Kt_rg=Kt_g                           # Kt_rg reduced global stiffness matrix

        #-------------Data from Bounday_Load_Conditions File--------------#
        #--------------------------------BCS------------------------------#
        #BCS:: An Array which is used to delete Global Stiffness matrix rows and coloumns
        #      for solving the equations
        Kt_rg=np.delete(Kt_rg,BCS,axis=0) #axis=0 is row
        Kt_rg=np.delete(Kt_rg,BCS,axis=1)
        reduced_G_global=G_global
        reduced_G_global=np.delete(reduced_G_global,BCS,axis=0)
        dU_g=np.matmul(np.linalg.inv(Kt_rg),-reduced_G_global)

        #-------For Newton Raphson Scheme covergence Criterion--------------#
        dU_g_convergence=dU_g
        u_g_convergence=np.delete(u_g,BCS,axis=0)
        dU_g_insert=dU_g
        #----dU_g_insert matrix is inserted with Boundary and Load conditions nodal values as 0------#
        for val in BCS:
            dU_g_insert=np.insert(dU_g_insert,val,0,axis=0)
        u_g=u_g+dU_g_insert

        if (np.linalg.norm(reduced_G_global,np.inf)<0.005*np.linalg.norm(F_g_int,np.inf) or np.linalg.norm(dU_g_convergence,np.inf)<0.005*np.linalg.norm(u_g_convergence,np.inf) or Newton > 5) :     
            break
        else:
            Newton=Newton+1
    #------------------End of "DO Newton_Raphson_method"----------------------------#
    if Newton>5:
        print("Convergence criterion is not met, at load (peercentage):",tau*100)
        break
    else:
        U_g_0 = u_g
print('No.of Newton raphson iterations done:',Newton)



            #---------------------------------------------------------------------#
            #                               Plotting                              #
            #---------------------------------------------------------------------#

#--------------------------Displacement Calculations for plotting----------------------------#

#Seperating out only mechanical displacements for plotting
Mechanical_Displacements=U_g_0[0:ncp*nudof]

#Reshape the mechanical displacements to match with shape of control points
Mechanical_Displacements = Mechanical_Displacements.reshape((ncpeta,ncpxi,2))

#Adding displacements to control points P_new = P + U  
#Pw is along with corresponding weights
Pw_new_temp = Pw[:,:,0:2] + Mechanical_Displacements

Pw_new = np.zeros_like(Pw)
Pw_new[:,:,0:2] = Pw_new_temp
Pw_new[:,:,3] = Pw[:,:,3]

# Transpose control points array is required for "NURBS_Surface_Point" Function 
Pw_new_T = Pw_new.transpose((1,0,2)) 
#-----------------------End of Displacement Calculations for plotting------------------------#


#------------------------Looping over NURBS surface points----------------------------#
x_positions=np.linspace(XI[0],(XI[-1]),ncpxi)
y_positions=np.linspace(ETA[0],(ETA[-1]),ncpeta)
surface=np.zeros((np.size(x_positions),np.size(y_positions),3))
for i,xi in enumerate(x_positions):
    for j,eta in enumerate(y_positions):
        # Here send input as displaced positions of control points matrix i.e P_new = P + XI
        S_r = NURBS_Surface_Point(n,p,XI,m,q,ETA,Pw_new_T,xi,eta)
        surface[i,j]=S_r
# Surface generated will be a transpose of what is desired since 
# NURBS_Surface_Point algorithm is modelled in such a way so a transpose is required.
surface = surface.transpose((1,0,2)) 
X=surface[:,:,0]
Y=surface[:,:,1]
#---------------------End of Looping over NURBS surface points--------------------------#


#------------------U1 Plotting----------------------#
fig,ax=plt.subplots(figsize = (7,5))
level = np.linspace(np.amin(Mechanical_Displacements[:,:,0]),np.amax(Mechanical_Displacements[:,:,0]),13)  
plt.contourf(X, Y, Mechanical_Displacements[:,:,0],levels=level,cmap='jet')
cbar = plt.colorbar(format='%.3e')
cbar.ax.tick_params(labelsize=15) 
plt.scatter(Pw_new[:,:,0],Pw_new[:,:,1])
ax.set(xlabel = 'x [$mm$]', ylabel = 'y [$mm$]', title= 'U1 (mm)')
ax.set_xlabel('x [$mm$]',fontsize = 15)
ax.set_ylabel('y [$mm$]', fontsize = 15)
ax.set_title('U1 (mm)',fontsize= 15)
output_file_path = Path("Results", "U1.png")
output_file_path.parent.mkdir(exist_ok=True)
fig.savefig(output_file_path)
#---------------End of U1 Plotting------------------#

#------------------U2 Plotting----------------------#
fig,ax=plt.subplots(figsize = (7,5))
level = np.linspace(np.amin(Mechanical_Displacements[:,:,1]),np.amax(Mechanical_Displacements[:,:,1]),13)
plt.contourf(X, Y, Mechanical_Displacements[:,:,1],levels=level,cmap='jet')
cbar = plt.colorbar(format='%.3e')
cbar.ax.tick_params(labelsize=15) 
plt.scatter(Pw_new[:,:,0],Pw_new[:,:,1])
ax.set_xlabel('x [$mm$]',fontsize = 15)
ax.set_ylabel('y [$mm$]', fontsize = 15)
ax.set_title('U2 (mm)',fontsize= 15)
output_file_path = Path("Results", "U2.png")
output_file_path.parent.mkdir(exist_ok=True)
fig.savefig(output_file_path)
#---------------End of U2 Plotting------------------#

#------------------Calculations for RF Plotting----------------------#
Mechanical_RF = F_g_int[0:ncp*nudof]
Mechanical_RF = Mechanical_RF.reshape((ncpeta,ncpxi,2))
RF1=Mechanical_RF[:,:,0]
RF2=Mechanical_RF[:,:,1]
#---------------End of Calculations for RF Plotting------------------#


#------------------RF1 Plotting----------------------#
fig,ax=plt.subplots(figsize = (7,5))
level = np.linspace(np.amin(RF1),np.amax(RF1),13)
plt.contourf(X, Y, RF1,levels=level,cmap='jet')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15) 
plt.scatter(Pw_new[:,:,0],Pw_new[:,:,1])
ax.set_xlabel('x [$mm$]',fontsize = 15)
ax.set_ylabel('y [$mm$]', fontsize = 15)
ax.set_title('RF1 (N)',fontsize= 15)
output_file_path = Path("Results", "RF1.png")
output_file_path.parent.mkdir(exist_ok=True)
fig.savefig(output_file_path)
#---------------------------------------------------#

#------------------RF2 Plotting----------------------#
fig,ax=plt.subplots(figsize = (7,5))
level = np.linspace(np.amin(RF2),np.amax(RF2),13)
plt.contourf(X, Y, RF2,levels=level,cmap='jet')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15) 
plt.scatter(Pw_new[:,:,0],Pw_new[:,:,1])
ax.set_xlabel('x [$mm$]',fontsize = 15)
ax.set_ylabel('y [$mm$]', fontsize = 15)
ax.set_title('RF2 (N)',fontsize= 15)
output_file_path = Path("Results", "RF2.png")
output_file_path.parent.mkdir(exist_ok=True)
fig.savefig(output_file_path)
#---------------------------------------------------#

#-----------------Calculations for Plotting Electrical Potential-----------------#
EPOT=U_g_0[nudof*ncp:(nudof*ncp+nedof*ncp)]
EPOT=EPOT.reshape((ncpeta,ncpxi))
#---------------------------------------------------------------------------------#

#-----------------Plotting Electrical Potential -----------------------------#
fig,ax=plt.subplots(figsize = (7,5))
level = np.linspace(np.amin(EPOT),np.amax(EPOT),13)
plt.contourf(X, Y, EPOT,levels=level,cmap='jet')
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15) 
plt.scatter(Pw_new[:,:,0],Pw_new[:,:,1])
ax.set_xlabel('x [$mm$]',fontsize = 15)
ax.set_ylabel('y [$mm$]', fontsize = 15)
ax.set_title('EPOT (Volts)',fontsize= 15)
output_file_path = Path("Results", "EPOT.png")
output_file_path.parent.mkdir(exist_ok=True)
fig.savefig(output_file_path)
#----------------------------------------------------------------------------#



#--------------Calculations for Plotting Electrical Reaction Forces---------------#
Electrical_RF = F_g_int[nudof*ncp:(nudof*ncp+nedof*ncp)]
Electrical_RF = -Electrical_RF.reshape((ncpeta,ncpxi))
#---------------------------------------------------------------------------------#

#-----------------Plotting Electrical Reaction Forces-------------------------#
fig,ax=plt.subplots(figsize = (7,5))
level = np.linspace(np.amin(Electrical_RF),np.amax(Electrical_RF),13)
plt.contourf(X, Y, Electrical_RF,levels=level,cmap='jet')
cbar = plt.colorbar(format='%.3e')
cbar.ax.tick_params(labelsize=15) 
plt.scatter(Pw_new[:,:,0],Pw_new[:,:,1])
ax.set_xlabel('x [$mm$]',fontsize = 15)
ax.set_ylabel('y [$mm$]', fontsize = 15)
ax.set_title('Reactive Electrical Charge (mC)',fontsize= 15)
output_file_path = Path("Results", "RCHG.png")
output_file_path.parent.mkdir(exist_ok=True)
fig.savefig(output_file_path)
#-----------------------------------------------------------------------------#

Main_end = timeit.timeit()
print('Time taken for Complete analysis to run is:',Main_end - Main_start)
sys.stdout.close()
sys.stdout=stdoutOrigin
