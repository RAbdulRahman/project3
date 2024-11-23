#Project 3
import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt


####### PART/QUESTION 1 #############



# to solve eqn 8 and 9 with B.C. m(r=0)=0 and rho (r =0) = rho c

#eqn 8: d rho/dr = -m rho/gamma(x)r^2
#eqn 9: dm/dr = r^2 rho , x = rho^1/3

def deriv1(r,initstate): 
    # defining derivative for part 1 using equations 8 and 9
    rho = initstate[0]
    m = initstate[1]
    if rho >= 0:
        x = (rho)**(1/3)
        gamma = x**2/(3*(1+x**2)**0.5)
        d_rho = -m*rho/(gamma*r**2)
        d_mass = (r**2) * rho
        return d_rho,d_mass


def stop_integ(r,state):
    # stopping when rho less than zero
    return state[0]

stop_integ.terminal = True
stop_integ.direction = -1


####   ANSWER 1 COMMENTS  #########

# The tiny value of r chosen here is 0.00000001. This would translate to about
# 0.000001 % of R0. Since this value is taken to be 7,720/ue km, the starting radius
# would be ~7/ue cm ~ 3.5 cm when ue =2.
# This small of a radius is negligibly different from starting at 0 for Earth sized bodies. 

#########
r_span = np.linspace(0.00000001,100,1000000)

#Testing Function 
initials = [2500000,0] # rho C, mass =0

# TEST PLOT
# result = scint.solve_ivp(deriv1,(0.00000001,100),initials,t_eval =r_span,events=stop_integ)

# plt.plot(result.t,result.y[0])
# plt.show()
# plt.plot(result.t,result.y[1])
# plt.show()
rhocs = [0.1,1,10,150,2000,7000,15000,80000,100000,150000,2500000]
i = 1
for rhoc in rhocs:
    initials = [rhoc,0]
    result = scint.solve_ivp(deriv1,(0.00000001,100),initials,t_eval =r_span,events=stop_integ)
    plt.subplot(3,4,i)
    plt.plot(result.t,result.y[0])
    plt.plot(result.t,result.y[1])
    plt.xlabel('Radius')
    plt.ylabel('Parameter')
    titlestrng = 'Rho C = ' + str(rhoc)
    plt.title(titlestrng)
    i += 1
    print('Max Mass (dimensionless)~',result.y[1][-1])
    
plt.legend(['Dimensionless Density','Dimensionless Mass'])
plt.suptitle('DImensionless Mass and Density as Funtions of Radius')
plt.tight_layout()
plt.show()



#####    PART/ANSWER 2 ########


# If ue chosen as 2 then M0 ~ 1.4175 e30 kg
# And R0 ~ 3.86 e6 m
# rho0 = 1.948 e9 kg/m^3

# Not using astropy since instructions specify 'MAY find ... useful'

MSun = 1.988400 * 10**30 #kg, source: https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
REarth = 6378000 #m, equatorial radius, source: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
m0 = 1.4175* 10**30 #kg
r0 = 3.86 * 10**6 #m
rho0 = 1.948 * 10**9 #kg/m^3
r_span = np.linspace(0.00000001,10,10000000)
for rhoc in rhocs: # Copying previous function but now with modified output for physicsal units
    initials = [rhoc*rho0,0]
    result = scint.solve_ivp(deriv1,(0.00000001,10),initials,t_eval =r_span,events=stop_integ) 
    plt.plot(result.t*r0,result.y[0]*rho0)
    plt.xlabel('Radius (m)')
    plt.ylabel('Density kg/m^3')
    plt.show()

    plt.plot(result.t*r0,result.y[1]*m0)
    plt.xlabel('Radius (m)')
    plt.ylabel('Mass (kg)')
    plt.show()

    
    print('Max Mass ()~',result.y[1][-1]*m0)




