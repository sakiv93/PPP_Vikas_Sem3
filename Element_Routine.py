#----------------------------------Displacement Driven------------------------------------------#
#--------------------------------------27th March-----------------------------------------------#
import numpy as np
from Geometry import *
from Material_Routine import *
from preprocessing_functions import *

def Jacobian12(xi,eta,elXI,elETA):
    """
    Input: 
        Input Parametric Co-ordinates (xi,eta) of Gauss Points
        ***Gauss Co-ordinates are defined in master space***
        Knot span range (elXI,elETA) in xi and eta direction
    Process: 
        Inputs functions involving partial derivatives of NURBS Basis Functions
        w.r.t to parametric co-ordinates and Calculates J1 matrix
    Return: 
        The function returns J1 and determinant of J1 and J2 matrix
    """
    #------------Calculating J1-------------------------#
    # J1 matrix for mapping from paramteric space to physical space.

    dxi_dximas = 0.5*(elXI[1]-elXI[0])          # Refer Equation 25 from Documentation
    deta_detamas = 0.5*(elETA[1]-elETA[0])      # Refer Equation 26 from Documentation 
    J2det = dxi_dximas*deta_detamas

    #------------Calculating J2-------------------------#
    # J2 matrix for mapping from master space to parametric space.

    #--------Evaluation of NURBS Surface basis functions derivatives---------#
    Aders,wders = SurfaceDerivsAlgAuv(n,p,XI,m,q,ETA,P,W,xi,eta,d)

    #-----------------------------------------------------------------------------------------------
    dRx_dxi     = (Aders[1][0][0]*wders[0][0] - wders[1][0]*Aders[0][0][0])/(wders[0][0]*wders[0][0])
    dRy_dxi     = (Aders[1][0][1]*wders[0][0] - wders[1][0]*Aders[0][0][1])/(wders[0][0]*wders[0][0])
    dRx_deta    = (Aders[0][1][0]*wders[0][0] - wders[0][1]*Aders[0][0][0])/(wders[0][0]*wders[0][0])
    dRy_deta    = (Aders[0][1][1]*wders[0][0] - wders[0][1]*Aders[0][0][1])/(wders[0][0]*wders[0][0])
    J1 = np.array([[dRx_dxi[0],dRx_deta[0]],
                [dRy_dxi[0],dRy_deta[0]]])

    J1det = (dRx_dxi[0]*dRy_deta[0])-(dRx_deta[0]*dRy_dxi[0])
    J1inv = np.linalg.inv(J1)
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
    xi_span = FindSpan(n,p,xi,XI)
    eta_span = FindSpan(m,q,eta,ETA)
    # Indices of non vanishing basis functions
    NVxi = np.arange(xi_span-p,xi_span+1,1)
    NVeta = np.arange(eta_span-q,eta_span+1,1)
    Aders,wders = SurfaceDerivsAlgAuv(n,p,XI,m,q,ETA,P,W,xi,eta,d)
    Denom = wders[0][0]
    Denom_du = wders[1][0]
    Denom_dv = wders[0][1]

    dR_dx = np.zeros((necpxi,necpeta))
    dR_dy = np.zeros((necpxi,necpeta))
    #---------------Loop over Non vanishing NURBS Basis Functions-------------#
    for ii, ii_value in enumerate(NVxi):
        for jj, jj_value in enumerate(NVeta):
            BFxi = BasisFuns(xi_span,xi,p,XI)
            BFeta = BasisFuns(eta_span,eta,q,ETA)
            Num = BFxi[ii]*BFeta[jj]**W[ii][jj]
            DBFxi = DersBasisFuns(xi_span,xi,p,d,XI)
            DBFeta = DersBasisFuns(eta_span,eta,q,d,ETA)
            Num_dxi = DBFxi[1][ii]*BFeta[jj]*W[ii][jj]
            Num_deta = BFxi[ii]*DBFeta[1][jj]*W[ii][jj]
            dR_dxi = Num_dxi/Denom - Denom_du*Num/(Denom*Denom)
            dR_deta = Num_deta/Denom - Denom_dv*Num/(Denom*Denom) 

            dR_dxi_dR_deta = np.array([dR_dxi,dR_deta])
            dR_dx_dR_dy = np.matmul(np.linalg.inv(J1),dR_dxi_dR_deta)
            dR_dx[ii][jj] = dR_dx_dR_dy[0]
            dR_dy[ii][jj] = dR_dx_dR_dy[1]
    #-------------Have to multiply dR_dx / dR_dy matrix with jacobian matrix-------------------#

    #---------Flatten (Convert 2D to 1D array) DervsNURBS Function------------#
    #dR_dx(0,0)....dR_dx(0,1),dR_dx(1,0),.....dR_dx(1,4).....dR_dx(4,0),..........dR_dx(4,4)
    fdR_dx = (np.transpose(dR_dx)).flatten()
    fdR_dy = (np.transpose(dR_dy)).flatten()
    #-------------------------Bu Matrix --------------------------#
    for i2 in range(necp):
        j1= 2*i2
        j2= 2*i2+1
        Bu[0,j1] = fdR_dx[i2]
        Bu[1,j2] = fdR_dy[i2]
        Bu[2,j1] = fdR_dy[i2]
        Bu[2,j2] = fdR_dx[i2]

        Be[0,i2] = fdR_dx[i2]
        Be[1,i2] = fdR_dy[i2]
    return Bu,Be

def elementRoutine(U_e,elXI,elETA,T_m):
    """
    Input: 
        Input (U_e) matrix which contain DOF values 
        Knot span range (elXI,elETA) in xi and eta direction
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

    Kt_e = np.zeros(((nudof+nedof)*necp,(nudof+nedof)*necp))            # Element stiffness matrix
    K_MM=np.zeros((nudof*necp,nudof*necp))                              # Refer Equation 56 from Documentation
    K_ME=np.zeros((nudof*necp,nedof*necp))                              # Refer Equation 57 from Documentation
    K_EM=np.zeros((nedof*necp,nudof*necp))                              # Refer Equation 58 from Documentation
    K_EE=np.zeros((nedof*necp,nedof*necp))                              # Refer Equation 59 from Documentation
    Fu_int_e=np.zeros((nudof*necp,1))                                   # Internal force matrix (Mechanical)
    Fe_int_e=np.zeros((nedof*necp,1))                                   # Internal force matrix (Electrical)
    F_int_e=np.zeros(((nudof+nedof)*necp,1))                            # Internal force matrix
    sigma_ig=np.zeros((np.shape(GPs_Ws)[0],3,1))                        # Stress matrix to store values at each gauss point
    epsilon_ig=np.zeros((np.shape(GPs_Ws)[0],3,1))                      # Strain matrix to store values at each gauss point
    electric_field_ig = np.zeros((np.shape(GPs_Ws)[0],2,1))             # Electrical field matrix to store values at each gauss point
    Electrical_Displacement_ig = np.zeros((np.shape(GPs_Ws)[0],2,1))    # Electrical displacement matrix to store values at each gauss point

    #-----------Looping over gauss point-----------#

    for j in range(np.shape(GPs_Ws)[0]):

        gp = GPs_Ws[j,0:2]      # For fetching respective gauss points
        wg = GPs_Ws[j,2]        # For fetching respective gauss point weights
        ximas = gp[0]       # Gauss points in Master space
        etamas = gp[1]      # Gauss points in Master space
        xi = 0.5*((elXI[1]-elXI[0])*ximas + (elXI[1]+elXI[0]))          # Gauss points in Parametric space
        eta = 0.5*((elETA[1]-elETA[0])*etamas + (elETA[1]+elETA[0]))    # Gauss points in Parametric space

        #------Calling Jacobian12 Function to get J1 and J2 determinants-----------#

        J1,J1det,J2det = Jacobian12(xi,eta,elXI,elETA)

        #-----------------Calling B_matrix function for Bu and Be Matrix--------------------#
        #-------Bu::(B matrix for mechanical case) Be::(B matrix for electrical case) ------#

        Bumatrix,Bematrix = B_matrix(xi,eta,J1)

        # U_u contains mechanical displcement from control point 1 to total control points in the element
        # U_phi contains electric potential from control point 1 to total control points in the element
        U_u=U_e[0:necp*nudof]
        U_phi = U_e[necp*nudof:]

        epsilon = np.matmul(Bumatrix,U_u)               # Strain
        electric_field = -np.matmul(Bematrix,U_phi)     # Electric_field
        epsilon_ig[j]=epsilon                   #Storing value of Strain at each gauss point
        electric_field_ig[j] = electric_field   #Storing value of electric_field at each gauss point

        #-----------------------Calling Matrial Routine----------------------#

        C, e, k, sigma, Electrical_Displacement = materialRoutine(epsilon,electric_field, T_m)
        #-----Storing value of Stress at each gauss point
        sigma_ig[j] = sigma 
        #-----Storing value of Electrical Displacement at each gauss point                                   
        Electrical_Displacement_ig[j] = Electrical_Displacement
        #-------------------------Local Stiffness matrix Ke-------------------#
        CBu=np.matmul(C,Bumatrix)
        BuCBu = np.matmul(np.transpose(Bumatrix),CBu)

        eBe=np.matmul(np.transpose(e),Bematrix)
        BueBe=np.matmul(np.transpose(Bumatrix),eBe)

        eBu=np.matmul(e,Bumatrix)
        BeeBu = np.matmul(np.transpose(Bematrix),eBu)

        kBe=np.matmul(k,Bematrix)
        BekBe = np.matmul(np.transpose(Bematrix),kBe)

        #------------------Numerical Integration-----------------------#

        K_MM = K_MM + BuCBu*J1det*J2det*wg*Thick
        K_ME = K_ME + BueBe*J1det*J2det*wg*Thick
        K_EM = K_EM + BeeBu*J1det*J2det*wg*Thick
        K_EE = K_EE - BekBe*J1det*J2det*wg*Thick

        #Arranging to Kt_e matrix (Local element Stiffness matrix)
        #    Kt_e = [K_MM K_ME]         # Refer Equation 55 from Documentation
        #           [K_EM K_EE]
        Kt_e[0:necp*nudof,0:necp*nudof]                                     = K_MM
        Kt_e[0:necp*nudof,necp*nudof:necp*(nudof+nedof)]                    = K_ME
        Kt_e[necp*nudof:necp*(nudof+nedof),0:necp*nudof]                    = K_EM
        Kt_e[necp*nudof:necp*(nudof+nedof),necp*nudof:necp*(nudof+nedof)]   = K_EE

        #------------Internal Force calculations--------------------#

        #-------Internal Force calculation for Mechanical case
        Fu_int_e= Fu_int_e+np.matmul(np.transpose(Bumatrix),sigma)*J1det*J2det*wg*Thick
        #-------Internal Force calculation for Electrical case
        Fe_int_e= Fe_int_e+np.matmul(np.transpose(Bematrix),Electrical_Displacement)*J1det*J2det*wg*Thick

        # ------Arranging to local Force_internal matrix 
        F_int_e[0:nudof*necp] = Fu_int_e                              # 0 1 2 3 .... 2*ncep-2, 2*ncep-1
        F_int_e[nudof*necp:(nudof*necp+nedof*necp)] = Fe_int_e        # 2*ncep, 2*ncep+1 .......  2*ncep+(ncep-1)
        #-------Initiating Force External as zero array because of displacement driven algorithm  
        F_ext_e = np.zeros_like(F_int_e)    
    return Kt_e, F_int_e, F_ext_e, sigma_ig, Electrical_Displacement_ig,epsilon_ig,electric_field_ig