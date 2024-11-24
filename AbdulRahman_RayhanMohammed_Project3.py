#Project 3
import numpy as np
import astropy.units as u
import scipy.integrate as scint
import matplotlib.pyplot as plt


####### PART/QUESTION 1 ###############################
#######################################################


# to solve eqn 8 and 9 with B.C. m(r=0)=0 and rho (r =0) = rho c

#eqn 8: d rho/dr = -m rho/gamma(x)r^2
#eqn 9: dm/dr = r^2 rho , x = rho^1/3

def deriv1(r,initstate): 
    # defining derivative for part 1 using equations 8 and 9
    '''Provide initial state as a list of density and mass'''
    rho = initstate[0]
    m = initstate[1]
    if rho >= 0:
        x = (rho)**(1/3)
        gamma = x**2/(3*(1+(x**2))**0.5)
        d_rho = -m*rho/(gamma*(r**2))
        d_mass = (r**2) * rho
        return d_rho,d_mass


def stop_integ(r,state):
    # stopping when rho less than zero
    return state[0]

stop_integ.terminal = True
stop_integ.direction = -1


####   ANSWER 1 COMMENTS  #########

# The tiny starting value of r chosen here is 0.00000001. This would translate to about
# 0.000001 % of R0. Since this value is taken to be 7,720/ue km, the starting radius
# would be ~7/ue cm ~ 3.5 cm when ue =2.
# This small of a radius is negligibly different from starting at 0 for Earth sized bodies. 

#########
r_span = np.linspace(0.00000001,100,1000000)

#Testing Function 
initials = [100,0] # rho C, mass =0

# Example plot
result = scint.solve_ivp(deriv1,(0.00000001,100),initials,t_eval =r_span,events=stop_integ)

plt.plot(result.t,result.y[0])
plt.plot(result.t,result.y[1])
plt.xlabel('Radius (multiple of R0)')
plt.ylabel('Mass or Desnity (multiple of M0 and ρ0 respectively)')
plt.title('Example Dimensionless Plot, ρ_c = 100')
plt.legend(['Density','Mass'])
plt.show()

# Ten rho c and their plot
rhocs = [0.1,0.5,1,2.5,5,10,25,100,500,2000,2500000]
i = 1
for rhoc in rhocs: # Solving over list of chosen rho_cs
    initials = [rhoc,0]
    result = scint.solve_ivp(deriv1,(0.00000001,100),initials,t_eval =r_span,events=stop_integ)
    if i == 9:
        i = 10 # to make sure legend doesnt cover graph
    plt.subplot(3,4,i) #plotting subplots on one figure
    plt.plot(result.t,result.y[0])
    plt.plot(result.t,result.y[1])
    plt.xlabel('Radius')
    plt.ylabel('Parameter')
    plt.title(f'ρ_c = {rhoc}') 
    i += 1 # chnaging position of subplot using i
    print('Mass (dimensionless)~',result.y[1][-1])
    
plt.legend(['Dimensionless Density','Dimensionless Mass'],loc='lower right')
plt.suptitle('Dimensionless Mass and Density as Funtions of Radius')
plt.tight_layout()
plt.show()



#######  PART/ANSWER 2 ########################################
###############################################################

# Constants, conversion done while plotting

r_sun = 695700000 #m , source: https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
m_sun = 1.988400 * 10**30 #kg
rho_sun = 1408 # mean desnity, kg/m3

# given values
rho0 = 1.948 * 10**9 #kg/m^3
m0 = 1.4175* 10**30 #kg
r0 = 1.544 * 10**7 /4 #m



r_span = np.linspace(0.00000001,10,10000000)

#seperate figures for physical mass and density, took out density figure for submission
#fig1 = plt.figure(1)
fig2 = plt.figure(2)

i=1

for rhoc in rhocs: # Copying previous function but now with modified output for physicsal units
    initials = [rhoc,0]
    result = scint.solve_ivp(deriv1,(0.00000001,10),initials,t_eval =r_span,events=stop_integ) 
    # plt.figure(1)
    # #plt.subplot(3,4,i)                              # Plotted rho to check if it made sense 
    # plt.plot(result.t*r0/r_sun,result.y[0]*rho0/rho_sun)
    # plt.xlabel('Radius (M⊙)')
    # plt.ylabel('Density (ρ⊙)')
    # plt.title(f'ρ_c = {rhoc}') 
    
    plt.figure(2) #plotting as instructed with mass on x axis and radius on y
    #plt.subplot(3,4,i)
    plt.plot(result.y[1]*m0/m_sun,result.t*r0/r_sun)
    plt.ylabel('Radius (R⊙)')
    plt.xlabel('Mass (M⊙)')
    
    i += 1

    print('Mass (M⊙)~',result.y[1][-1]*m0/m_sun, 'for', 'Radius (R⊙)~',result.t[-1]*r0/r_sun)

# plt.figure(1)
# plt.suptitle('Density vs Radius for Various Starting Densities')
# plt.tight_layout()


plt.figure(2) 
plt.suptitle('Radius vs Mass for Various Starting Densities')
plt.tight_layout()
plt.show()



########## ANSWER/PART 2 COMMENTS #########

# Using Kippenhahn & Weigert's Chandrasekhar limit definition, the Chandrasekhar limit is 5.836/4 * Msun =  1.459 M⊙.
# As can be seen in the (maximum) masses printed, our estimate would be about 1.433 M⊙ which is 
# reasonably close to the Kippenhahn & Weigert citation with a percentage difference of around 1.8%. 
# This difference could be attributed to integration errors and improvements in the accuracy of sun's mass since 1990. 





######## PART/ANSWER 3 ###################################################
##########################################################################

rho_cs = [1,2,5]
for rho_c in rho_cs:
    initials = [rho_c,0]
    result = scint.solve_ivp(deriv1,(0.00000001,100),initials,t_eval =r_span,events=stop_integ) # Default method is RK45
    result2 = scint.solve_ivp(deriv1,(0.00000001,100),initials,t_eval =r_span,events=stop_integ,method = 'DOP853') # DOP853 should be more precise than standard RK 45
    # plt.plot(result.y[1]*m0/m_sun,result.t*r0/r_sun,)
    # plt.xlabel('Mass (M⊙)')
    # plt.ylabel('Radius (R⊙)')
    print('Mass (M⊙) using RK45: ',round(result.y[1][-1]*m0/m_sun,4),'Mass (M⊙) using DOP853: ',round(result2.y[1][-1]*m0/m_sun,4))
    print('Radius (R⊙) using RK45: ',round(result.t[-1]*r0,4),'Radius (R⊙) using DOP853: ',round(result2.t[-1],4))
    print()








### PART 4 ################################################################
###########################################################################


# Plotting Mass vs Radius, rho Cs chosen to fit csv 
rho_cs = [1,1.5,2,3,5]
for rho_c in rho_cs:
    initials = [rho_c,0]
    result = scint.solve_ivp(deriv1,(0.00000001,100),initials,t_eval =r_span,events=stop_integ)
    plt.plot(result.y[1]*m0/m_sun,result.t*r0/r_sun,)
    plt.xlabel('Mass (M⊙)')
    plt.ylabel('Radius (R⊙)')
   

# lists for csv data
x = []
x_err = []
y = []
y_err = []

with open('wd_mass_radius.csv', 'r') as file:
    lines = file.readlines()[1:]  # Ignore the first line (header)
    
    for line in lines: # filling lists to plot csv data
        values = line.strip().split(',')
        x.append(float(values[0]))        
        x_err.append(float(values[1]))    
        y.append(float(values[2]))        
        y_err.append(float(values[3]))    

# Plot data with error bars
plt.errorbar(x, y, xerr=x_err, yerr=y_err, fmt='o', capsize=2, label='Data')


plt.legend(['ρ_c = 1','ρ_c = 1.5','ρ_c = 2','ρ_c = 3','ρ_c = 5','Given Data'],loc = 'upper right')
plt.title('Plot for Chosen Values of ρC ')
plt.show()





