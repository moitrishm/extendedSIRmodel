import numpy as np
import random 
import sys
import time
#from scipy.integrate import odeint
#import matplotlib.pyplot as plt
from math import *
#import scipy.io as sio
#from scipy import signal

t_i = time.time()
# Total population, N.
N1 = 100000
N2 = 5000000
# Initial number of infected and recovered individuals, I0 and R0.
I10=100
R10=0
I20 =1000
# Everyone else, S0, is susceptible to infection initially.
S10 = N1 - I10 - R10
S20 = N2 - I20
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta_hv=0.03
alpha=0.0011
gamma=0.01
mu_h = 0.02
                                                                                                                     
beta_vh =0.05
mu_0 = 0.125
mu_v = 0.125
a=float(sys.argv[1])

# A grid of time points (in days)
Tf0 = 5
delta = 365*Tf0
t0 = np.linspace(0, delta-1, delta)

R0  = np.zeros(np.int(delta))
#mu= np.zeros(np.int(delta))

for i in range(np.int(delta)):
   R0[i] = (N2*np.exp(mu_0*(1+a*np.cos(2*np.pi*t0[i]/365))- mu_v)*beta_hv*beta_vh)/(N1*(gamma+mu_h)*mu_0 * (1 +  a * np.cos(2*np.pi*t0[i]/365)) )
   #mu[i] = mu_v * (1 +  a * np.cos(2*np.pi*t[i]/365))
            

# The SIR model differential equations.
#def deriv(y, t0, N2, beta_hv, alpha, gamma, mu_h, beta_vh, mu_v):
 #   S1, I1, R1, S2, I2 = y
    ####
  #  dS1dt = mu_h*N1 - beta_hv * I2 * S1 / (S1+I1+R1) + alpha * R1 - mu_h * S1
   # dI1dt = beta_hv * I2 * S1 / (S1+I1+R1) - gamma * I1 - mu_h * I1 
    #dR1dt = gamma * I1 - alpha * R1 - mu_h * R1                          
  #  dS2dt = mu_0 *(1+a*np.cos(2*np.pi*t0/365))*(S2+I2) - beta_vh * I1 * S2 / (S1+I1+R1) - mu_v* S2
   # dI2dt = beta_vh * I1 * S2 / (S1+I1+R1) - mu_v * I2      
    #return dS1dt, dI1dt, dR1dt, dS2dt, dI2dt 

# Initial conditions vector
#y0 = S10, I10, R10, S20, I20

# Integrate the SIR equations over the time grid, t.
#sol = odeint(deriv, y0, t0, args=(N2, beta_hv, alpha, gamma, mu_h, beta_vh, mu_v))
#Sh, Ih, Rh, Sv, Iv = sol.T #transpose


#t1=t0/365


# Final time (no. of iterations of simulation)
Tf=int(9*1e7)
Tf1 =  int(Tf/10)
#R0  = np.zeros(Tf)
T  = np.zeros(Tf1)
S1 = np.zeros(Tf1)
I1 = np.zeros(Tf1)
R1 = np.zeros(Tf1)
S2 = np.zeros(Tf1)
I2= np.zeros(Tf1)

t = 0.0


T[0] = t 
S1[0] = S10
I1[0] = I10
R1[0] = R10
S2[0] = S20
I2[0] = I20


tS1, tI1, tR1, tS2, tI2 = S10, I10, R10, S20, I20
  
count = 0
while((tI1 > -0.01 ) and (tI2 > -0.01) and (count < Tf-1) and t < 100.0 ):   #t is number of days for which simulation will take place
    #if(t%365 < )
    k1 = (mu_0)*(1+a*np.cos(2*np.pi*t/365))*(tS2 +tI2)
    k2 = mu_v* tS2
    k3 = mu_v * tI2
    k4 = gamma * tI1
    k5 =  mu_h * tI1 
    k6 = mu_h * tR1
    k7 = alpha * tR1
    k8 = beta_vh * tI1 * tS2 / (tS1+tI1+tR1)
    k9 = mu_h * N1
    k10 = mu_h * tS1
    k11 = beta_hv * tI2 * tS1 / (tS1+tI1+tR1)
        

    K = k1 + k2 + k3 + k4 + k5 +k6 + k7 + k8 + k9 +k10 + k11
    #R0[count] = (mu_v)*(1+a*np.cos(2*np.pi*t/365))
    delta = (1.0/K) * np.log(1.0/random.random())
    t += delta
    r = random.random() * K
    if (r < k1)  :
        tS2 += 1
    elif (r < k1+k2) :
        tS2 -= 1
    elif (r < k1+k2+k3) :
        tI2 -= 1
    elif (r < k1+k2+k3+k4) :
        tI1 -= 1
        tR1 += 1
    elif (r < k1+k2+k3+k4+k5) :
        tI1 -= 1
    elif (r < k1+k2+k3+k4+k5+k6) :
        tR1 -= 1
    elif (r < k1+k2+k3+k4+k5+k6+k7) :
        tR1 -= 1
        tS1 += 1
    elif (r < k1+k2+k3+k4+k5+k6+k7+k8) :
        tS2 -= 1
        tI2 += 1
    elif (r < k1+k2+k3+k4+k5+k6+k7+k8+k9) :
        tS1 += 1
    elif (r < k1+k2+k3+k4+k5+k6+k7+k8+k9+k10) :
        tS1 -= 1
    elif (r < k1+k2+k3+k4+k5+k6+k7+k8+k9+k10+k11) :
        tS1 -= 1
        tI1 += 1
           
    
    
    #print(  round(k1,2),  round(k2,2), round(k3,2),round(k4,2),round(k5,2),round(k6,2),round(k7,2),round(k8,2),round(k9,2),round(k10,2),round(k11,2), round(r,2), tS1,tI1,tR1,tS2,tI2 )
    count += 1    
    if ((count%10)==0):
      ind=int(count/10)
    
      T[ind]  = t
      S1[ind] = tS1
      I1[ind] = tI1
      R1[ind] = tR1
      S2[ind] = tS2
      I2[ind] = tI2
    

print(  round(tI1,2),  round(tI2,2), count, ind)
#R0 = R0[:count]
T = T[:ind]
S1 = S1[:ind]
I1 = I1[:ind]
R1 = R1[:ind]
S2 = S2[:ind]
I2 = I2[:ind]


np.savez_compressed('output_{}'.format(a) , a=T, b=I1, c=I2, d=R1, e=S2) 
t_f = time.time() - t_i
print( t_f)

