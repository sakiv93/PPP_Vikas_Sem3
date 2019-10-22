import numpy as np 
import math

def N(i,k,u):
    if t(i) <= u <= t(i+1):
        N=1
    else:
        N=0
    return N

#Parameters
control_points=9
degree=3
u1=control_points-degree+2
u1=control_points-degree+2
u=np.arange(0,u1,1)

u1=np.linspace(i,i+1,10)
for i in range(control_points+1):
    for u in u1:
        N[i,k] = ((u-t[i])*(N(i,k-1)/(t[i+k-1]-t[i])) \
            + ((t[i+k]-u)*(N(i+1,k-1)/(t[i+k]-t[i+1]))

#Extracting knot values
control_points=5
degree=3
t = np.zeros(control_points+degree+1)
for i in range (control_points+degree+1):
    if i<degree:
        t[i]=0
    elif degree<=i<=control_points:
        t[i]=i-degree+1
    else:
        t[i]=control_points-degree+2

