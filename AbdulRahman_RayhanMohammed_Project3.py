#Project 3
import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt


# Question 1

# to solve eqn 8 and 9 with B.C. m(r=0)=0 and rho (r =0) = rho c

#eqn 8: d rho/dr = -m rho/gamma(x)r^2
#eqn 9: dm/dr = r^2 rho , x = rho^1/3

def deriv1(r,initstate): 
    rho = initstate[0]
    m = initstate[1]
    x = (rho)**(1/3)
    gamma = x**2/(3*(1+x**2)**0.5)

    d_rho = -m*rho/(gamma*r**2)
    d_mass = r**2 * rho

    return d_rho,d_mass


def stop_integ(r,state):
    return state[0]

stop_integ.terminal = True
stop_integ.direction = -1

r_span = np.linspace(0.00001,10,10000)
initials = [0.1,0] # rho C, mass =0

test = scint.solve_ivp(deriv1,(0.00001,10),initials,t_eval =r_span,events=stop_integ)

plt.plot(test.t,test.y[0])
plt.show()
plt.plot(test.t,test.y[1])
plt.show()

