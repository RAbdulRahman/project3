#Project 3
import numpy as np
import scipy.integrate as scint


# Question 1

# to solve eqn 8 and 9 with B.C. m(r=0)=0 and rho (r =0) = rho c

#eqn 8: d rho/dr = -m rho/gamma(x)r^2
#eqn 9: dm/dr = r^2 rho , x = rho^1/3

def deriv1(r,initstate): 
    rho = initstate[0]
    m = initstate[1]
    x = rho**(1/3)
    gamma = x**2/(3*(1+x**2)**0.5)

    d_rho = -m*rho/(gamma*r**2)
    d_mass = r**2 * rho



