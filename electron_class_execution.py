from electron_class import Electron
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import random
import math
e1 = Electron()

'''
#executing pressure 
print(len(e1.p(x = np.linspace(50000, 0, 5000))))
'''
#-------------------------------
'''
#executing Z
x = np.linspace(50000, 0, 5000)
e1.p(x)
print(len(e1.Z(x,373,200)))
'''
#-------------------------------
'''
#executing refractive index
x = np.linspace(50000, 0, 5000)
e1.p(x)
e1.Z(x,373,200)
print(e1.n(0.3,x,15,0,450))
'''
#-------------------------------
'''
#executing plot n vs h
x = np.linspace(50000, 0, 5000)
plt.plot(e1.n(0.3,x,15,0,450),x)
#plt.plot(x, e1.n(0.3,x,15,0,450))
plt.ylabel('height(cm)')
plt.xlabel('n')
plt.show()
e1.p(x)
e1.Z(x,373,200)
e1.n(0.3,x,15,0,450)
#e1.plot_n_h()
'''
#-------------------------------
'''
#executing gammas for all gamma0 
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
e1.bremsstrahlung(x)
e1.total_momentum_100(x)
print(e1.gamma1_new)
'''
#-------------------------------
'''
#solution for bremsstrahlung
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
e1.bremsstrahlung(x)
e1.total_momentum_100(x)
#print(e1.gamma1_new)
#print(e1.brem1)
print(e1.sol1)
'''
#---------------------------
'''
#executing beta for gamma
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
#print(e1.beta_electron())
e1.beta_electron()
plt.plot(x, e1.sol1_beta1, 'b', label='$\gamma_{0}$=1000')
plt.show()
'''
#-------------------------------
'''
#executing gammas for minimum and maximum wavelength
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
e1.beta_electron()
print(e1.eta_gamma(x,0.3,0.7))
'''
#-------------------------------right
'''
#executing plot for gamma vs x(Ionization+Bremstrahlung)
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
#e1.bremsstrahlung(x)
#e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
#e1.plot_gamma_h(x)
plt.plot(x/10**2, e1.sol1, 'b', label=r'$\gamma_{0}$=100', linewidth=2.0)
#plt.plot(x, e1.brem1, 'b--', label=r'$\gamma_{0}$=100')---------------------------------------only Bremstrahlung which is not need according to the script
plt.plot(x/10**2, e1.sol2, 'g', label=r'$\gamma_{0}$=1500', linewidth=2.0)
#plt.plot(x, e1.brem2, 'g--', label=r'$\gamma_{0}$=1500')--------------------------------------only Bremstrahlung '' '' ''  ''    ''     ''     ''     ''
plt.plot(x/10**2, e1.sol3, 'r', label=r'$\gamma_{0}=10^4$', linewidth=2.0)
#plt.plot(x, e1.brem3, 'r--', label=r'$\gamma_{0}=10^4$')-------------------------------------only Bremstrahlung ''       ''            ''           ''
plt.axhline(y=160, color='c', label=r'$E_{c}=160$', linewidth=2.0)
plt.plot(x/10**2, e1.gamma1, 'k--', label=r'$\gamma_{min}$($\lambda$=300)', linewidth=2.0) # --------------lorentz factor corresponding to wavelength 300nm
plt.plot(x/10**2, e1.gamma2, 'm--', label=r'$\gamma_{min}$($\lambda$=700)', linewidth=2.0) # ------------lorentz factor corresponding to wavelength 700nm
plt.ylim(30, 10000)
#plt.xlim(49730,50000)#-----------------------------------------------limits to get a proper view of lines if we add energy loss of electrons due to Bremstrahlung and Ionization
for var in (e1.sol1, e1.sol2, e1.sol3):
    plt.annotate('%0.2f' % var.max(), xy=(1, var.max()), xytext=(8, 0), 
                 xycoords=('axes fraction', 'data'), textcoords='offset points')
plt.xlabel('Height(m)')
plt.ylabel(r'Lorentz factor($\gamma$)')
plt.yscale('log')
plt.rcParams.update({'font.size': 22})
plt.grid()
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.14), ncol=3, fancybox=True, shadow=True)
plt.show()
'''
#------------------------------------right
'''
#-----------------------------------wrong

#total gamma from Ionization and bremsstrahlung
#x = np.linspace(50000, 0, 5000)
#e1.gamma_lorentz(x)
#e1.bremsstrahlung(x)
#e1.eta_gamma(x,0.3,0.7)
#e1.total_momentum_100(x)
#e1.total_momentum_500(x)
#e1.total_momentum_1000(x)
#e1.beta_electron()
#print(e1.sol1)
#print(e1.brem1)
#print(e1.gamma1_new)
#plt.plot(x, e1.gamma1_new, 'b', label=r'$\gamma_{0}$=100')
#plt.plot(x, e1.gamma2_new, 'g', label=r'$\gamma_{0}$=500')
#plt.plot(x, e1.gamma3_new, 'r', label=r'$\gamma_{0}=10^3$')
#plt.plot(x, e1.gamma1, 'k--', label=r'$\gamma_{min}$($\lambda$=300)') # --------------lorentz factor corresponding to wavelength 300nm
#plt.plot(x, e1.gamma2, 'm--', label=r'$\gamma_{min}$($\lambda$=700)') # ------------lorentz factor corresponding to wavelength 700nm
#plt.ylim(30, 10000)
#plt.xlim(49730,50000)#-----------------------------------------------limits to get a proper view of lines if we add energy loss of electrons due to Bremstrahlung and Ionization
#plt.xlabel('height(cm)')
#plt.ylabel(r'$\gamma$')
#plt.yscale('log')
#plt.grid()
#plt.rcParams.update({'font.size': 16})
#plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.14), ncol=3, fancybox=True, shadow=True)
#plt.show()
#---------------------------------------------------wrong
'''



#-------------------------------
'''
#executing refractive index for blue and red
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
print(e1.nref_b_r())
'''
#-------------------------------
'''
#executing plot for angle of emission vs x
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
#e1.bremsstrahlung(x)
e1.eta_gamma(x,0.3,0.7)
#e1.total_momentum_100(x)
#e1.total_momentum_500(x)
#e1.total_momentum_1000(x)
e1.beta_electron()
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_1500()
e1.conversion_10000()
e1.conversion_red()
e1.beta_n()
#e1.plot_angle_h(x)
plt.plot(x/10**2,e1.theta3, 'b', label='$\lambda_{0}$=300($\gamma_{0}$=100)',linewidth=2.0)
plt.plot(x/10**2,e1.theta4, 'b--', label='$\lambda_{0}$=700($\gamma_{0}$=100)',linewidth=2.0)
plt.plot(x/10**2,e1.theta3_1500, 'g', label='$\lambda_{0}$=300($\gamma_{0}$=1500)',linewidth=2.0)
plt.plot(x/10**2,e1.theta4_1500, 'g--', label='$\lambda_{0}$=700($\gamma_{0}$=1500)',linewidth=2.0)
#plt.plot(x/10**2,e1.theta3_10000, 'r', label='$\lambda_{0}$=300($\gamma_{0}=10^{4}$)',linewidth=2.0)
#plt.plot(x/10**2,e1.theta4_10000, 'r--', label='$\lambda_{0}$=700',linewidth=2.0)
#plt.xlim(-0.01,50000)
plt.ylim(0,1.6)
plt.rcParams.update({'font.size': 22})
plt.xlabel('Height(m)')
plt.ylabel(r'Angle of emission:${\Theta_{c}}(^{0})$')
leg = plt.legend()
leg_lines = leg.get_lines()
leg_texts = leg.get_texts()
#plt.setp(leg_lines, linewidth=6)
#plt.setp(leg_texts, fontsize='x-large')
#plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.14), ncol=3, fancybox=True, shadow=True)
plt.legend(loc='lower left')
plt.grid()
plt.show()
'''

#-------------------------------
'''

#executing for radius of electron rings
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
print(e1.radius(x))
'''
#------------------------------
'''
#executing plot for radius vs x
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
#e1.bremstrahlung(x)
e1.eta_gamma(x,0.3,0.7)
#e1.total_momentum_100(x)
#e1.total_momentum_500(x)
#e1.total_momentum_1000(x)
e1.beta_electron()
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_1500()
e1.conversion_10000()
e1.conversion_red()
e1.beta_n()
e1.radius(x)
#e1.plot_radius_h(x)
plt.plot(x/10**2,e1.radius_b/10**2, 'b', label='$\lambda_{0}$=300($\gamma_{0}$=100)',linewidth=2.0)
plt.plot(x/10**2,e1.radius_r/10**2, 'b--', label='$\lambda_{0}$=700($\gamma_{0}$=100)',linewidth=2.0)
plt.plot(x/10**2,e1.radius_b_1500/10**2, 'g', label='$\lambda_{0}$=300($\gamma_{0}$=1500)',linewidth=2.0)
plt.plot(x/10**2,e1.radius_r_1500/10**2, 'g--', label='$\lambda_{0}$=700($\gamma_{0}$=1500)',linewidth=2.0)
#plt.plot(x/10**2,e1.radius_b_10000/10**2, 'r', label='$\lambda_{0}$=300($\gamma_{0}=10^{4}$)',linewidth=2.0)
#plt.plot(x/10**2,e1.radius_r_10000/10**2, 'r--', label='$\lambda_{0}$=700',linewidth=2.0)
plt.ylim(0,15)
plt.rcParams.update({'font.size': 22})
plt.xlabel('Height(m)')
plt.ylabel('Radius(m)')
#plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
plt.legend(loc='upper left')
plt.grid()
plt.show()
'''
#-----------------------------------
'''
#executing for number of potons per path length(dndx_num)
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
e1.radius(x)
print(e1.dndx_num(x))
'''
#----------------------------------------

#executing plot for number of photons emitted per path length(dN/dx vs x)
'''
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
#e1.bremstrahlung(x)
e1.eta_gamma(x,0.3,0.7)
#e1.total_momentum_100(x)
#e1.total_momentum_500(x)
#e1.total_momentum_1000(x)
e1.beta_electron()
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_1500()
e1.conversion_10000()
e1.dndx_num(x)
e1.dndx_num_1500(x)
e1.dndx_num_10000(x)
plt.plot(x/10**2,e1.Result3, 'b', label='$\gamma_{0}$=100',linewidth=2.0)
plt.plot(x/10**2,e1.Result3_1500, 'g', label='$\gamma_{0}$=1500',linewidth=2.0)
plt.plot(x/10**2,e1.Result3_10000, 'r', label='$\gamma_{0}$=10000',linewidth=2.0)
#plt.xlabel('Height(m)')
plt.xlabel('x(m)')
plt.ylabel(r'${dN/dx}(m^{-1})$')
#plt.xlim(0,50000)
plt.rcParams.update({'font.size': 22})
plt.ylim(0,50)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fancybox=True, shadow=True)
plt.grid()
plt.show()
'''
#--------------------------------------

'''
#executing total number of photons from dN/dx integrating over x(N=19992359.2906)

x = np.linspace(500.00, 0, 5000)
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
e1.radius(x)
e1.dndx_num(x)
print(e1.N_dndx(x))
'''
#---------------------------------
'''
#calculating angle of emission using a range of wavelength l(0.3-0.7 micrometer):
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
e1.radius(x)
print(e1.theta(x))
#e1.theta(x)
#plt.plot(x,e1.theta0)
#plt.show()
'''
#-----------------------------------
'''
#calculating slope of the curve theta_ch vs x taking the derivative
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
e1.theta(x)
print(np.absolute(e1.slope_dtheta_dx(x)))
#print(len(e1.slope_dtheta_dx(x)))
'''
'''
#plot angle vs slope:
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
e1.theta(x)
e1.slope_dtheta_dx(x)
plt.plot(e1.theta0, e1.func1)
plt.show()
'''
'''
#calculating the dN/domega
x = np.linspace(50000, 0, 5000)
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
e1.theta(x)
e1.slope_dtheta_dx(x)
print(len(e1.dndomega(x)))
'''
'''
#plotting dN/dOmega vs theta_ch
x = np.linspace(50000, 0, 5000)#in cm
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
e1.theta(x)
e1.slope_dtheta_dx(x)
e1.dndomega(x)
#print(e1.theta0*(180/np.pi))
plt.plot(e1.theta0*(180/np.pi),e1.Result4)
plt.xlabel(r'${\Theta}$(degree)')
plt.ylabel(r'${dN/{d\Omega}}$')
#plt.xlim(0.0795720656263863152,1.27)#limits of theta values
plt.legend()
plt.grid()
plt.show()
'''
'''
#plotting dN/dOmega vs height:
x = np.linspace(50000, 0, 5000)#in cm
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
e1.theta(x)
e1.slope_dtheta_dx(x)
e1.dndomega(x)
plt.plot(x,e1.Result4, 'b', label='$\gamma_{0}$=1e4')
plt.xlabel('height(cm)')
plt.ylabel(r'${dN/{d\Omega}}$')
plt.legend()
plt.grid()
plt.show()
'''

#printing no per pixel
'''
x = np.linspace(50000, 0, 5000)#in cm
phi = np.linspace(0, 360, 5000)#l is the \phi coordinate, which range from 1 to 360(in degrees). 
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
e1.theta(x)
e1.slope_dtheta_dx(x)
e1.dndomega(x)
e1.u_v_conversion(phi)
print(e1.Result4)
'''


#converting dN/d\Omega vs \theta_ch to u-v image

x = np.linspace(50000, 0, 5000)#in cm
phi = np.random.uniform(0, 360, 5000)#l is the \phi coordinate, which range from 1 to 360(in degrees). 
phi1 = phi*math.pi/180#--------------------------\phi in radian to make the calculation easy(since angle u and v used radian)
e1.gamma_lorentz(x)
#e1.bremstrahlung(x)
e1.eta_gamma(x,0.3,0.7)
#e1.total_momentum_100(x)
#e1.total_momentum_500(x)
#e1.total_momentum_1000(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_1500()
e1.conversion_10000()
e1.conversion_red()
e1.beta_n()
e1.theta(x)
e1.slope_dtheta_dx(x)
e1.dndomega(x)
e1.u_v_conversion(phi1)
#plt.scatter(np.absolute(e1.u),np.absolute(e1.v))
#plt.hist2d(e1.u, e1.v, bins=(50, 50), norm=LogNorm(), normed=True, weights = e1.Result4, range=[[-0.1,0.1],[-0.1,0.1]])
plt.hist2d(e1.u, e1.v, bins=(60, 60), norm=LogNorm(), normed=True, weights = e1.Result4, range=[[-0.1,0.1],[-0.1,0.1]])
plt.grid()
cb = plt.colorbar()
cb.set_label('Number of cherenkov photons/pixel')
plt.rcParams.update({'font.size': 22})
plt.xlabel('u')
plt.ylabel('v')
plt.show()
#print(e1.v)



#theta from u-v(hessII) , Height of emission(both in histogram)
'''
x = np.linspace(50000, 0, 5000)#-----------------height of interest
phi = np.linspace(0, 360, 5000)#-----------------\phi coordinate, which range from 1 to 360(in degrees). 
phi1 = phi*math.pi/180#--------------------------\phi in radian to make the calculation easy(since angle u and v used radian)
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
e1.theta(x)
e1.u_v_conversion(phi1)
#print(len(e1.value))
#plt.hist(x, bins=50,range = None,log = True)#, normed = True)
plt.hist(e1.theta0*(180/np.pi), bins=50,range = None)#,log = True)
plt.xlabel('angle of emission(degrees)')
#plt.xlabel('emission height(cm)')
plt.show()
'''

'''
x = np.linspace(50000, 0, 5000)#in cm
e1.gamma_lorentz(x)
e1.beta_electron()
e1.eta_gamma(x,0.3,0.7)
e1.nref_b_r()
e1.conversion_blue()
e1.conversion_red()
e1.beta_n()
e1.theta(x)
e1.slope_dtheta_dx(x)
e1.dndomega(x)
e1.u_v_conversion()
e1.circle()
plt.contour(e1.X,e1.Y,e1.F,[0])
plt.show()
'''
