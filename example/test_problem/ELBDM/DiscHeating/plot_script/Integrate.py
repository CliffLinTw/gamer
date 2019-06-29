from scipy import integrate
import numpy as np
import math
import matplotlib.pyplot as plt

TotM = 6.743
G =  4.43738E-6
B  = 0.00875    ### A exp(-r/B)
Rmax = 0.0175
div = 10
A  = TotM/(2*np.pi*B*(B-(Rmax+B)*np.exp(-Rmax/B)))

m_range = 30  ### Series expansion
L_range = 30  ### Series expansion

Data = np.zeros([div+1,6])  ### 0:radius, 1: potential by scipy, 2: error by scipy, 3: potential by series
                                 ### 4: acc by series, 5: vel by series

###########  Scipy

for i in range(div):
    r = Rmax/div*(i+1)

    def I(r_prime, theta):
        density = A*np.exp(-r_prime/B)
        distance = np.sqrt(r**2+r_prime**2-2*r*r_prime*np.cos(theta))
        return G * density*r_prime/distance

    Data[i+1,0] = r                                         #### Radius
    answer = integrate.nquad(I, [[0, Rmax],[0, 2*np.pi]])
    Data[i+1,1] = answer[0]                                 #### Potential by Scipy
    Data[i+1,2] = answer[1]                                 #### Error of Integration
    print('r =',Data[i+1,0],'potential =',Data[i+1,1])

##########  functions of Series Expansion:

leading = 2*np.pi*G*A   ## Leading term

def factorial(num):
    if num == 0:
        return 1
    else:
        fac = 1
        for i in range(num):
            fac *= (i+1)
        return fac

def L_coe(L_input):    ###
    if L_input == 0:
        return 2*np.pi
    else:
        L_coeff = 2*np.pi
        for Label in range(L_input):
            L_coeff *= ((2*Label+1)/(2*Label+2))**2
        return L_coeff

def Left_part(L,m,r):  ### for acc
    if L-m == 1:
        return r**m*((m+1)/(L+m+2))
    else:
        return r**m*((m+1)/(L+m+2)+(m+1)/(L-m-1))

def Right_part(L,m,r):  ### for acc
    if L-m == 1:
        return (r**(L-1))*(-1+L*np.log(Rmax/r))
    else:
        return -L*(r**(L-1))*(Rmax**(-L+m+1))/(L-m-1)

def Left_part_phi(L,m,r):  ### for phi
    if L-m == 1:
        return r**(m+1)/(L+m+2)
    else:
        return r**(m+1)*(1/(L+m+2)+1/(L-m-1))

def Right_part_phi(L,m,r):  ### for phi
    if L-m == 1:
        return (r**L)*(np.log(Rmax/r))
    else:
        return -(r**L)*(Rmax**(-L+m+1))/(L-m-1)

########### Series expansion

for i in range(div):
    r = (Rmax/(div-1)*i)
    valueA = 0
    valuephi = 0
    for m in range(m_range):
        m_part = 1/factorial(m)*(-1/B)**m
        for L in range(L_range):
            if r == 0 :
                valueA += 0
            else:
                L_part = L_coe(L)
                LP = Left_part(L,m,r)
                RP = Right_part(L,m,r)
                LPP = Left_part_phi(L,m,r)
                RPP = Right_part_phi(L,m,r)
                valueA   += leading*m_part*L_part*(LP+RP)
                valuephi += leading*m_part*L_part*(LPP+RPP)

    Data[i,4] = valueA
    Data[i,3] = valuephi
    if valueA > 0 :
        Data[i,5] = 0
    else:
        Data[i,5] = 1.0*math.sqrt((abs(valueA*r)))

####### plot and output

plt.scatter(Data[:,0],Data[:,1]) ### phi by scipy
plt.scatter(Data[:,0],Data[:,3]) ### phi by series
#plt.scatter(Data[:,0],Data[:,4]) ### acc by series
#plt.scatter(Data[:,0],Data[:,5]) ### vel by series

#plt.ylim(-0.2E-34, 1.5E-34)
plt.show()

#np.savetxt("profile.txt", potential)
