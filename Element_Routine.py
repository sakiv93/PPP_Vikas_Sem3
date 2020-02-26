#----------------------------------Displacement Driven------------------------------------------#
#------------------------------------------26 Feb-----------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
np.set_printoptions(threshold=np.inf)
from Geometry import *
from Material_Routine import *
from preprocessing_functions import *

def Jacobian12(xi,eta,elU,elV):
    """
    Input: 
        Input Parametric Co-ordinates (xi,eta) of Gauss Points
        ***Gauss Co-ordinates are defined in master space***
        Knot span range (elU,elV) in xi and eta direction
    Process: 
        Inputs functions involving partial derivatives of NURBS Basis Functions
        w.r.t to parametric co-ordinates and Calculates J1 matrix
    Return: 
        The function returns Determinant of J1 and J2 matrix
    """
    #------------Calculating J1-------------------------#

    dxi_dximas = 0.5*(elU[1]-elU[0])
    deta_detamas = 0.5*(elV[1]-elV[0]) 
    J2det = dxi_dximas*deta_detamas

    #------------Calculating J2-------------------------#

    #--------Evaluation of NURBS Surface basis functions derivatives---------#
    Aders,wders = SurfaceDerivsAlgAuv(n,p,U,m,q,V,P,W,xi,eta,d)
    #print('Aders:',Aders)
    #print('wders:',wders)

    #-----------------------------------------------------------------------------------------------
    #---------------------- Without using RatSurfaceDerivs----------------------------#
    dRx_dxi     = (Aders[1][0][0]*wders[0][0] - wders[1][0]*Aders[0][0][0])/(wders[0][0]*wders[0][0])
    dRy_dxi     = (Aders[1][0][1]*wders[0][0] - wders[1][0]*Aders[0][0][1])/(wders[0][0]*wders[0][0])
    dRx_deta    = (Aders[0][1][0]*wders[0][0] - wders[0][1]*Aders[0][0][0])/(wders[0][0]*wders[0][0])
    dRy_deta     = (Aders[0][1][1]*wders[0][0] - wders[0][1]*Aders[0][0][1])/(wders[0][0]*wders[0][0])
    J1 = np.array([[dRx_dxi[0],dRx_deta[0]],
                [dRy_dxi[0],dRy_deta[0]]])
    print('J1',J1)


    #-----------------------------------------------------------------------------------------------
    # dRx = RatSurfaceDerivs(Aders,wders,d)
    # #print('dRx',dRx)
    # dRx_dxi  =  dRx[1][0][0]
    # dRx_deta =  dRx[0][1][0]
    # dRy_dxi  =  dRx[1][0][1]
    # dRy_deta =  dRx[0][1][1]
    # print('dRx',dRx[1][0][0],dRx[0][1][0],dRx[1][0][1],dRx[0][1][1])
    # J1 = np.array([[dRx_dxi,dRx_deta],
    #             [dRy_dxi,dRy_deta]])
    # print('J1 Matrix : ',J1)
    #-----------------------------------------------------------------------------------------------


    J1det = (dRx_dxi[0]*dRy_deta[0])-(dRx_deta[0]*dRy_dxi[0])
    J1inv = np.linalg.inv(J1)
    #print('J1inv : ',J1inv)
    return J1,J1det,J2det

def B_matrix(xi,eta,J1):
    """
    Input: 
        Input Parametric Co-ordinates (xi,eta) of Gauss Points
        ***Gauss Co-ordinates are defined in master space***
    Process: 
        Calls function to find derivatives of NURBS Basis functions
        and Store the values in Bu matrix (B matrix related to mechanical case)
        and Be matrix (B matrix related to electrical case)
    Return: 
        The function returns Bu matrix and Be matrix
    """

    Bu=np.zeros((3,nudof*necp))
    Be=np.zeros((2,necp))
    #---------------NURBS Basis Functions Derivatives wrt x and y-------------#
    uspan = FindSpan(n,p,xi,U)
    print('uspan',uspan)
    vspan = FindSpan(m,q,eta,V)
    print('vspan',vspan)
    # Indices of non vanishing basis functions
    NVu = np.arange(uspan-p,uspan+1,1) #Thoroughly checked. Can crosscheck again
    print('NVu',NVu)
    NVv = np.arange(vspan-q,vspan+1,1)
    print('NVv',NVv)
    Aders,wders = SurfaceDerivsAlgAuv(n,p,U,m,q,V,P,W,xi,eta,d)
    Denom = wders[0][0]
    Denom_du = wders[1][0]
    Denom_dv = wders[0][1]
    #***********Change this back again************************************************
    # dR_du = np.zeros((ncpxi,ncpeta))
    # dR_dv = np.zeros((ncpxi,ncpeta))
    #***********************************************************

    dR_dx = np.zeros((necpxi,necpeta))
    dR_dy = np.zeros((necpxi,necpeta))
    #---------------Loop over Non vanishing NURBS Basis Functions-------------#
    #***************Forgot to multiply with control points weights*************# Have to do it
    for ii, ii_value in enumerate(NVu):
        for jj, jj_value in enumerate(NVv):
            BFu = BasisFuns(uspan,xi,p,U)
            print('BasisFuncsU',BFu)
            BFv = BasisFuns(vspan,eta,q,V)
            print('BasisFuncsV',BFv)
            # Num = BFu[ii]*BFv[jj]
            Num = BFu[ii]*BFv[jj]**W[ii][jj]
            print('Num',Num)
            DBFu = DersBasisFuns(uspan,xi,p,d,U)
            DBFv = DersBasisFuns(vspan,eta,q,d,V)
            Num_du = DBFu[1][ii]*BFv[jj]*W[ii][jj]
            Num_dv = BFu[ii]*DBFv[1][jj]*W[ii][jj]
            #Num_du = DBFu[1][ii]*BFv[jj]
            #Num_dv = BFu[ii]*DBFv[1][jj]
            #***********Change this back again************************************************
            # dR_du[ii_value][jj_value] = Num_du/Denom - Denom_du*Num/(Denom*Denom)
            # dR_dv[ii_value][jj_value] = Num_dv/Denom - Denom_dv*Num/(Denom*Denom)
            #***********Change this back again************************************************
            dR_dxi = Num_du/Denom - Denom_du*Num/(Denom*Denom)
            dR_deta = Num_dv/Denom - Denom_dv*Num/(Denom*Denom) 

            dR_dxi_dR_deta = np.array([dR_dxi,dR_deta])
            print('dR_dxi_dR_deta',dR_dxi_dR_deta)
            dR_dx_dR_dy = np.matmul(np.linalg.inv(J1),dR_dxi_dR_deta)
            print('dR_dx_dR_dy',dR_dx_dR_dy)
            dR_dx[ii][jj] = dR_dx_dR_dy[0]
            dR_dy[ii][jj] = dR_dx_dR_dy[1]
            
            # dR_dxi[ii][jj] = Num_du/Denom - Denom_du*Num/(Denom*Denom)
            # dR_deta[ii][jj] = Num_dv/Denom - Denom_dv*Num/(Denom*Denom)
            print('dR_du',dR_dx[ii][jj])
            print('dR_dv',dR_dy[ii][jj])
    #-------------Have to multiply dR_du / dR_dv matrix with jacobian matrix-------------------#






    #---------Flatten (Convert 2D to 1D array) DervsNURBS Function------------#
    #dR_du(0,0)....dR_du(0,1),dR_du(1,0),.....dR_du(1,4).....dR_du(4,0),..........dR_du(4,4)
    #print(dR_du)
    #fdR_du = dR_du.flatten()
    fdR_dx = (np.transpose(dR_dx)).flatten() #************Cross check this******************#
    print('fdR_dx',fdR_dx)
    #print(fdR_du)
    #fdR_dv = dR_dv.flatten()
    fdR_dy = (np.transpose(dR_dy)).flatten()
    print('fdR_dv',fdR_dy)
    #print(fdR_dv)
    #-------------------------Bu Matrix --------------------------#
    for i2 in range(necp):
        print('necp',necp)
        j1= 2*i2
        j2= 2*i2+1
        Bu[0,j1] = fdR_dx[i2]
        #print(j1,j2)
        #print(Bu[0,j1])
        Bu[1,j2] = fdR_dy[i2]
        #print(Bu[1,j2])
        Bu[2,j1] = fdR_dy[i2]
        #print(Bu[2,j1])
        Bu[2,j2] = fdR_dx[i2]
        #print(Bu[2,j2])

        Be[0,i2] = fdR_dx[i2]
        Be[1,i2] = fdR_dy[i2]
    return Bu,Be

def elementRoutine(U_e,elU,elV,T_m):
    """
    Input: 
        Input (U_e) matrix which contain DOF values 
        Knot span range (elU,elV) in xi and eta direction
    Process: 
        The element routine takes in information about the element 
        and loops over gauss points and perform numerical integration
        to find the values of Internal and External force matrix,
        Elemental stiffness matrix.
    Return: 
        The function returns Elemental Stiffness matrix (Kt_e), 
        Internal and External elemental force matrix, 
        Stress,Strain,Electric field and Electric displacements at Gauss points
    """
    print('Elu,Elv:',elU,elV)
    #Bu=np.zeros((3,nudof*necp))
    #Be=np.zeros((2,necp))

    Kt_e = np.zeros(((nudof+nedof)*necp,(nudof+nedof)*necp))
    K_MM=np.zeros((nudof*necp,nudof*necp))
    K_ME=np.zeros((nudof*necp,nedof*necp))
    K_EM=np.zeros((nedof*necp,nudof*necp))
    K_EE=np.zeros((nedof*necp,nedof*necp))
    Fu_int_e=np.zeros((nudof*necp,1))
    Fe_int_e=np.zeros((nedof*necp,1))
    F_int_e=np.zeros(((nudof+nedof)*necp,1))
    sigma_ig=np.zeros((np.shape(GPs_Ws)[0],3,1))
    epsilon_ig=np.zeros((np.shape(GPs_Ws)[0],3,1))
    electric_field_ig = np.zeros((np.shape(GPs_Ws)[0],2,1))
    Electrical_Displacement_ig = np.zeros((np.shape(GPs_Ws)[0],2,1))
    #F_ext_e=np.zeros(((nudof+nedof)*necp,1)

    #-----------Looping over gauss point-----------#

    for j in range(np.shape(GPs_Ws)[0]):

        gp = GPs_Ws[j,0:2]
        wg = GPs_Ws[j,2]
        ximas = gp[0]  #Gauss points in Master space
        etamas = gp[1] #Gauss points in Master space
        xi = 0.5*((elU[1]-elU[0])*ximas + (elU[1]+elU[0]))     #Gauss points in Parametric space
        eta = 0.5*((elV[1]-elV[0])*etamas + (elV[1]+elV[0]))   #Gauss points in Parametric space

        #------Calculating J1 and J2 determinants-----------#

        J1,J1det,J2det = Jacobian12(xi,eta,elU,elV)
        print('Determinat of J1 : ',J1det)
        print('Determinat of J2 : ',J2det)
        print('Weight',wg)

        #-------------------------Bu and Be Matrix --------------------------#

        Bumatrix,Bematrix = B_matrix(xi,eta,J1)
        print('xi,eta',xi,eta)
        print('Bu:',Bumatrix)
        print('Be:',Bematrix)

        # U_u contains mechanical displcement from node 1 to 4, and U_phi contains electric potential from node 1 to 4
        U_u=U_e[0:necp*nudof]
        U_phi = U_e[necp*nudof:]
        print('U_u',U_u)
        print('U_phi',U_phi)

        epsilon = np.matmul(Bumatrix,U_u)               # Strain
        electric_field = -np.matmul(Bematrix,U_phi)     # Electric_field
        epsilon_ig[j]=epsilon                   #Storing value of Strain at each gauss point
        electric_field_ig[j] = electric_field   #Storing value of electric_field at each gauss point
        print('epsilon',epsilon)
        print('electric_field',electric_field)

        C, e, k, sigma, Electrical_Displacement = materialRoutine(epsilon,electric_field, T_m)
        print('Sigma,ED',sigma,Electrical_Displacement)
        sigma_ig[j] = sigma
        Electrical_Displacement_ig[j] = Electrical_Displacement
        #-------------------------Local Stiffness matrix Ke-------------------#
        CBu=np.matmul(C,Bumatrix)
        BuCBu = np.matmul(np.transpose(Bumatrix),CBu)

        eBe=np.matmul(np.transpose(e),Bematrix)
        BueBe=np.matmul(np.transpose(Bumatrix),eBe)

        eBu=np.matmul(e,Bumatrix)
        BeeBu = np.matmul(np.transpose(Bematrix),eBu)

        kBe=np.matmul(k,Bematrix)
        BekBe = -np.matmul(np.transpose(Bematrix),kBe)

        #------------------Numerical Integration-----------------------#

        K_MM = K_MM + BuCBu*J1det*J2det*wg*Thick
        K_ME = K_ME + BueBe*J1det*J2det*wg*Thick
        K_EM = K_EM + BeeBu*J1det*J2det*wg*Thick
        K_EE = K_EE + BekBe*J1det*J2det*wg*Thick
        #print('KEE:',K_EE)

        #Arranging to Kt_e matrix 
        Kt_e[0:necp*nudof,0:necp*nudof]                                 = K_MM
        Kt_e[0:necp*nudof,necp*nudof:necp*(nudof+nedof)]                = K_ME
        Kt_e[necp*nudof:necp*(nudof+nedof),0:necp*nudof]                 = K_EM
        Kt_e[necp*nudof:necp*(nudof+nedof),necp*nudof:necp*(nudof+nedof)] = K_EE
        print('K_MM',Kt_e[0:necp*nudof,0:necp*nudof])
        print('K_ME',Kt_e[0:necp*nudof,necp*nudof:necp*(nudof+nedof)])
        print('K_EM',Kt_e[necp*nudof:necp*(nudof+nedof),0:necp*nudof])
        print('K_EE', Kt_e[necp*nudof:necp*(nudof+nedof),necp*nudof:necp*(nudof+nedof)])

        Fu_int_e= Fu_int_e+np.matmul(np.transpose(Bumatrix),sigma)*J1det*J2det*wg*Thick
        Fe_int_e= Fe_int_e+np.matmul(np.transpose(Bematrix),Electrical_Displacement)*J1det*J2det*wg*Thick

        #Arranging to F_int matrix 
        F_int_e[0:nudof*necp] = Fu_int_e                              # 0 1 2 3 4 5 6 7
        F_int_e[nudof*necp:(nudof*necp+nedof*necp)] = Fe_int_e        # 8 9 10 11
        F_ext_e = np.zeros_like(F_int_e)    

    #$$$$$     Didnt return gauss points co-ordinates     $$$$$#    
    return Kt_e, F_int_e, F_ext_e, sigma_ig, Electrical_Displacement_ig,epsilon_ig,electric_field_ig


#----------------Test case------------------#

# U_e = np.zeros((necp*nudof,1))
# Kt_e, F_int_e, F_ext_e,sigma = elementRoutine(U_e, 1)

# print(Kt_e)
# print(F_int_e)
# print(F_ext_e)
# print(sigma)