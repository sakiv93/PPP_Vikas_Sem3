#----------------------------------Displacement Driven------------------------------------------#
#--------------------------------------27th March-----------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import special
from mpl_toolkits.mplot3d import Axes3D  
np.set_printoptions(threshold=np.inf)

#------------------------------------------------------------------------------------------------#
#---------------Knot vector---------------#
def KnotVector(degree,n_control_points):
    """
    Input: 
        Input Degree of the NURBS curve (degree) 
        and number of control points (n_control_points)
    Process: 
        Function calculates Knot vector array for given 
        Degree and no.of Control points of NURBS Curve
    Return: 
        The function returns knot vector array
    """
    knot_vector = np.zeros((n_control_points-1)+(degree+1)+1)
    for i in range ((n_control_points-1)+(degree+1)+1):
        if i<(degree+1):
            knot_vector[i]=0
        elif (degree+1)<=i<=(n_control_points-1):
            knot_vector[i]=i-(degree+1)+1
        else:
            knot_vector[i]=(n_control_points-1)-(degree+1)+2
    return(knot_vector)
#------------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------------#
#---------------Find Span---------------#
# Function to find knot span in Knot vector
# Example:  If Knot vector = [0,0,1,2,3,4,5] and if u = 1.5 
#           then knot span is 2 (i.e index of 1 in knot vector) since 1.5 lies between 1 and 2 
def FindSpan(n_inp,degree_inp,u_inp,knot_vector_inp):
    """
    Input: 
        Input [(no.control points-1) - (degree+1)] as (n_inp) 
        degree of NURBS curve (degree_inp),
        Parametric Co-ordinate(u_inp) on the curve
        knot vector of the curve (knot_vector_inp)
    Process: 
        The function performs linear search for the knot span in knot vector
    Return: 
        The function returns the span of the parametric co-ordinate in the knot vector
    """
    x=knot_vector_inp[degree_inp+1]                 # Second Unique knot
    y=knot_vector_inp[-(degree_inp+2)]              # Last but one Unique knot
    if (u_inp < x ):
        return degree_inp
    elif (u_inp>y):
        return (len(knot_vector_inp)-(degree_inp+2))
    else:
        for i,pos in enumerate(knot_vector_inp):
            if math.floor(u_inp) == pos:
                return (i)
#------------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------------#
#-----------------Basis Functions---------------#
# Adapted from alogorithm A2.2 in NURBS Book page no.70
#The output is stored in N[0],....N[p]
#These are non zero functions in a given span
#i is the span of the u value
#we get N[i-p].........,N[i]
# So if i is 3 then N[3] is stored in N[p] of the output matrix
def BasisFuns(i,u,p,U):
    """
    Input: 
        Input knot span (i) of parametric co-ordinate (u) on NURBS curve
        Degree of the curve (p) and knot vector of the NURBS curve (U)

    Process: 
        The function find Non-zero basis function and also avoid 
        the occurence of divide by zero
        The output is stored in N[0],....N[p]
        These are non zero functions in a given span
        i is the span of the u value
        we get N[i-p].........,N[i]
        So if i is 3 then N[3] is stored in N[p] of the output matrix
    Return: 
        The function returns non zero basis functions for the given parametric co-ordinate
    """
    N =np.zeros(p+1)
    N[0] = 1.0
    left =np.zeros(p+2)
    right =np.zeros(p+2)
    for j in range(1,p+1):
        left[j] = u - U[i+1-j]
        right[j] = U[i+j] - u
        saved= 0.0
        for r in range(0,j):
            temp = N[r]/(right[r+1]+left[j-r])
            N[r] = saved + right[r+1]*temp
            saved = left[j-r]*temp
        N[j] = saved
    return N

#-----------------Derivatives of Basis Functions---------------#
# Adapted from alogorithm A2.3 in NURBS Book page no.72
# Derivatives are stored in ders[k][j]
# ders[k][j] is the kth derivative of function N_(i-p+j,p)
# If k>p the derivatives are zero.
def DersBasisFuns(i,u,p,m,U):
    """
    Input: 
        Input knot span (i) of parametric co-ordinate (u) on NURBS curve
        Degree of the curve (p), derivatives upto and including (m) th 
        knot vector of the NURBS curve (U)

    Process: 
        The function finds Non-zero basis function and also their derivatives
        Derivatives are stored in ders[k][j]
        ders[k][j] is the kth derivative of function N_(i-p+j,p)
        If k>p the derivatives are zero.
    Return: 
        The function returns non zero basis functions and their derivatives
        upto and including mth derivative for the given parametric co-ordinate
    """
    #Inititation of dimentions for 2D matrices
    ndu=np.zeros((p+1,p+1))
    ders=np.zeros((p+1,p+1))
    a=np.zeros((2,p+1))
    left =np.zeros(p+2)
    right =np.zeros(p+2)
    
    ndu[0][0]=1.0
    for j in range(1,p+1):
        left[j] = u - U[i+1-j]
        right[j] = U[i+j] - u
        saved=0.0
        for r in range(j):
            #Lower triangle
            ndu[j][r] = right[r+1]+left[j-r]
            temp=ndu[r][j-1]/ndu[j][r]
            #Upper triangle
            ndu[r][j] = saved+(right[r+1]*temp)
            saved=left[j-r]*temp
        ndu[j][j] = saved
    for j in range (p+1): #Load the basis functions
        ders[0][j]=ndu[j][p]
    #This secion computes the derivatives
    for r in range(p+1):
        s1=0
        s2=1 #Alternative rows in array a
        a[0][0] = 1.0
        #Loop to compute kth derivative
        for k in range(1,m+1):
            d=0.0
            rk=r-k
            pk=p-k
            if(r>=k):
                a[s2][0]=a[s1][0]/ndu[pk+1][rk]
                d=a[s2][0]*ndu[rk][pk]
            if(rk>=-1):
                j1=1
            else:
                j1=-rk
            if(r-1<=pk):
                j2=k-1
            else:
                j2=p-r
            for j in range (j1,j2+1):
                a[s2][j] =(a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j]
                d += (a[s2][j]*ndu[rk+j][pk])
            if(r<=pk):
                a[s2][k]=-a[s1][k-1]/ndu[pk+1][r]
                d+=(a[s2][k]*ndu[r][pk])
            ders[k][r]=d
            #Switch rows
            j=s1
            s1=s2
            s2=j
            #Multiply through by the correct factors
    r=p
    for k in range(1,m+1):
        for j in range(p+1):
            ders[k][j] =ders[k][j]* r
        r =r* (p-k)
    return ders

#--------------------------Aders,Wders-----------------------------------#

#A modified version of Algorithm A3.6 from NURBS Book Page no. 111
#An extra W is given as input to function for convinience which contain weights of corresponding control points
def SurfaceDerivsAlgAuv(n,p,U,m,q,V,P,W,u,v,d):
    #print('P,W,AdersWders',P,W)
    SKL_A=np.zeros((20,20,3))
    SKL_W=np.zeros((20,20,1))
    temp_A=np.zeros((20,3))
    temp_W=np.zeros((20,1))
    du = min(d,p)
    for k in range(p+1,d+1):
        for l in range(d-k+1):
            SKL_A[k][l]=0.0
            SKL_W[k][l]=0.0
    dv = min(d,q)
    for l in range(q+1,d+1):
        for k in range(d-l+1):
            SKL_A[k][l]=0.0 
            SKL_W[k][l]=0.0 
    uspan=FindSpan(n,p,u,U)
    unders=DersBasisFuns(uspan,u,p,du,U)
    vspan=FindSpan(m,q,v,V)
    vnders=DersBasisFuns(vspan,v,q,dv,V)
    for k in range(du+1):
        for s in range(q+1):
            temp_A[s] = 0.0
            temp_W[s] = 0.0
            for r in range(p+1):
                temp_A[s]=temp_A[s] + unders[k][r]*P[uspan-p+r][vspan-q+s]*W[uspan-p+r][vspan-q+s]
                temp_W[s]=temp_W[s] + unders[k][r]*W[uspan-p+r][vspan-q+s]
            dd=min(d-k,dv)
            for l in range(dd+1):
                SKL_A[k][l] = 0.0
                SKL_W[k][l] = 0.0
                for s in range(q+1):
                    SKL_A[k][l] = SKL_A[k][l] + vnders[l][s]*temp_A[s]
                    SKL_W[k][l] = SKL_W[k][l] + vnders[l][s]*temp_W[s]
    return SKL_A,SKL_W

#----------------Point on NURBS Surface---------------------#
def NURBS_Surface_Point(n,p,U,m,q,V,Pw,u,v):
    """
    Input: 
        Input [(no.control points-1) - (degree+1)] in xi and eta direction as {n} and {m}
        Degree of the curve in xi and eta direction {p} and {q} respectively
        knot vector of the NURBS curve in xi and eta direction {U} and {V} respectively
        Control point matrix as {Pw} 
        parametric co-ordinates {u} and {v} on NURBS curve in xi and eta direction 
        respectively
    Process: 
        The function performs tensor product between NURBS basis functions
        in xi and eta directions.
    Return: 
        The function returns physical co-ordinates of the NURBS surface point from the 
        given parametric co-ordinates {u} and {v}
    """
    surface_point=np.zeros(3)
    uspan = FindSpan(n,p,u,U)
    Nu=BasisFuns(uspan,u,p,U)
    vspan = FindSpan(m,q,v,V)
    Nv=BasisFuns(vspan,v,q,V)
    S=np.zeros(4)
    for d in range(4):
        Sw=0.0
        for j in range(q+1):
            for i in range(p+1):
                Sw = Sw + Nu[i]*Nv[j]*Pw[uspan-p+i][vspan-q+j][d]
        S[d]=Sw   
    #surface_point=S[:-1]
    surface_point=S[:-1]/S[-1]
    return surface_point

#****** Algorithm from IGA Simplified paper ********#

#-----------------Inputs-----------------------#
# n-no.of control points along xi direction
# p-Degree of basis function along xi direction
# m-no.of control points along eta direction
# q-Degree  of basis function along eta direction
def ControlPointAssembly(n,p,m,q,ele_no):
    nel= (n-p)*(m-q) 
    ncp= n*m
    necp = (p+1)*(q+1)
    ControlPointAssemblyArray = np.zeros((necp,nel))

    A=0
    e=0
    for j in range(1,m+1):
        for i in range(1,n+1):
            A=A+1
            if (i>=(p+1) and j>=(q+1)):
                e=e+1
                for jloc in range(q+1):
                    for iloc in range(p+1):
                        B=A-jloc*n-iloc
                        b=jloc*(p+1)+iloc+1
                        ControlPointAssemblyArray[b-1,e-1]=B
    CP = np.flip(np.transpose(ControlPointAssemblyArray),axis=1)-1 #To generate data in proper order
    # -1 so that indices will start from '0'
    print('CP',CP)
    return  CP[ele_no]