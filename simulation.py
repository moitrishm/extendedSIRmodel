import numpy as np
import random 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from math import *
import scipy.io as sio
from scipy import signal
import time

t_i = time.time()
# Total population, N.
N1 = 5000
N2 = 50000
# Initial number of infected and recovered individuals, I0 and R0.
I10=1000
R10=10
I20 =1000
# Everyone else, S0, is susceptible to infection initially.
S10 = N1 - I10 - R10
S20 = N2 - I20
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
beta_hv=0.03
alpha=0.0011
gamma=0.01
mu_h = 0.02
                                                                                                                     
beta_vh =0.1
mu_0 = 0.1
mu_v = 0.1
a=0.75    
Tf=int(9*1e8)

#Tf0 = 1
#delta = 365*Tf0
t = 0.0


#t0 = np.linspace(0, delta-1, delta)

#R0  = np.zeros(np.int(delta))
#mu= np.zeros(np.int(delta))

#for i in range(np.int(delta)):
   #R0[i] = (N2*np.exp (mu_0*(1+(signal.square(2*np.pi*t/365, duty=a)))- mu_v)*beta_hv*beta_vh)/(N1*(gamma+mu_h)*mu_0 * (1 + (signal.square(2*np.pi*t/365, duty=a))) )
#   mu[i] = (mu_0/a)*(1+signal.square(2*np.pi*t0[i]/365, duty=a))

#t1=t0/365


#R0  = np.zeros(Tf)
T  = np.zeros(Tf)
S1 = np.zeros(Tf)
I1 = np.zeros(Tf)
R1 = np.zeros(Tf)
S2 = np.zeros(Tf)
I2= np.zeros(Tf)




T[0] = t 
S1[0] = S10
I1[0] = I10
R1[0] = R10
S2[0] = S20
I2[0] = I20


tS1, tI1, tR1, tS2, tI2 = S10, I10, R10, S20, I20
  
count = 0
while((tI1 > -0.01 ) and (tI2 > -0.01) and (count < Tf-1) and t < 3650.0 ):   #t is number of days for which simulation will take place
    #if(t%365 < )
    #k1 = (mu_0/a)*(1+(signal.square(2*np.pi*t/365, duty=a)))*(tS2 +tI2)
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
	
  
	
    if ((count%5)==0): 
    	T[count]  = t
    	S1[count] = tS1
    	I1[count] = tI1
    	R1[count] = tR1
    	S2[count] = tS2
    	I2[count] = tI2
    

print(  round(tI1,2),  round(tI2,2), count, T[count])
#R0 = R0[:count]
T = T[:count]
S1 = S1[:count]
I1 = I1[:count]
R1 = R1[:count]
S2 = S2[:count]
I2 = I2[:count]

np.savez_compressed('output_medium'.format(a) , a=T, b=I1, c=I2, d=R1, e=S2) 
t_f = time.time() - t_i
print( t_f)



