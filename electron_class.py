'''
class for the calculation of total number of photons emitted per pixel
'''
#importing modules required for the calculation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from scipy.integrate import odeint
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from scipy import integrate
from scipy.integrate import quad, romberg
import math
from numpy import diff
import random

#defining constants used in the calculation:all these are defined inside __init__

#x = np.linspace(50000, 0, 20000)------------- height in the atmosphere concerned in cm
#K = 0.307-----------------MeV g^-1 cm^2
#zsq = 1.-------------------elementary charge of electron
#at_no = 7.5 ------------Atomic number of absorbers
#Am = 15.-------------Average atomic mass of oxygen and nitrogen(only major components components are considered) 
#En = 0.5110----------------Electron mass * c^2 MeV
#Isq = np.square(10.*at_no*1e-12) --------------mean excitation energy in eV
#h0 = 840000. ------------------Mean height in cm
#rho0 = 1e-3  ---------------density on the surface(g/cm^3}
#psurf = 1.   --------------------- surface pressure in bar
#alpha0 = 0.0073  ---------------Fine structure constant
#Rad_len = 37.1  #-------------Radiation length for dry air in g/cm^2
#initial conditions: for y = \gamma to solve homogeneous differential equation:defined in the __init__ methods:
#y0 = 40.
#y1 = 120.
#y2 = 70.
#y3 = 80.
#y4 = 50.

class Electron:
	def __init__(self):     
		self.y1 = 100
		self.y2 = 1500
		self.y3 = 10000
		self.K = 0.307
		self.zsq = 1.
		self.at_no = 7.5
		self.Am = 15.
		self.En = 0.5110
		self.h0 = 840000
		self.rho0 = 1e-3
		self.psurf = 1
		self.alpha0 = 0.0073
		self.Rad_len = 31e3  #-------------Radiation length for electrons in dry air is 37.1 g/cm^2 converted to cm in correspondence with the height(cm)
		#self.theta_pixel_hess = 0.00066    #--------------------half angle pixel of hess phase II
#defining pressure as a function of height
	def p(self, x):
		return (self.psurf*np.exp(-x/self.h0))*10**5  #----------------------10**5 is added to convert bar to Pa

#method to calculate the compressibility--------------------------------------------
#Reference:https://github.com/polyanskiy/refractiveindex.info-scripts/blob/master/scripts/Ciddor%201996%20-%20air.py
    	def Z(self,x,T,xw): 
        	t=T-273.15
        	a0 = 1.58123e-6   
        	a1 = -2.9331e-8   
        	a2 = 1.1043e-10   
        	b0 = 5.707e-6     
        	b1 = -2.051e-8    
        	c0 = 1.9898e-4    
        	c1 = -2.376e-6    
        	d = 1.83e-11     
        	e = -0.765e-8    
        	return(1-(self.p(x)/T)*(a0+a1*t+a2*t**2+(b0+b1*t)*xw+(c0+c1*t)*xw**2) + (self.p(x)/T)**2*(d+e*xw**2))

#defining refractive index------------------------------------------------------------
#Reference:https://github.com/polyanskiy/refractiveindex.info-scripts/blob/master/scripts/Ciddor%201996%20-%20air.py
	def n(self,l,x,t,h,xc):
        	sigma = 1/l           
        	T= t + 273.15    
        	R = 8.314510      
        	k0 = 238.0185     
        	k1 = 5792105      
        	k2 = 57.362       
        	k3 = 167917       
        	w0 = 295.235      
        	w1 = 2.6422       
        	w2 = -0.032380    
        	w3 = 0.004028     
        	A = 1.2378847e-5  
        	B = -1.9121316e-2 
        	C = 33.93711047
        	D = -6.3431645e3  
        	alpha = 1.00062
        	beta = 3.14e-8       
        	gamma = 5.6e-7        

        	if(t>=0):   #------------------------------------saturation vapor pressure of water vapor in air at temperature T
            		svp = np.exp(A*T**2 + B*T + C + D/T) 
        	else:
            		svp = 10**(-2663.5/T+12.537)
    
        	f0 = alpha + beta*self.p(x) + gamma*t**2
        	xw = f0*h*svp/self.p(x)
        	nas = 1 + (k1/(k0-sigma**2)+k3/(k2-sigma**2))*1e-8
        	naxs = 1 + (nas-1) * (1+0.534e-6*(xc-450))
        	nws = 1 + 1.022*(w0+w1*sigma**2+w2*sigma**4+w3*sigma**6)*1e-8
        	Ma = 1e-3*(28.9635 + 12.011e-6*(xc-400)) 
        	Mw = 0.018015                            
        	Za = self.Z(x,288.15, 0)                
        	Zw = self.Z(x,293.15, 1)                  
        	rhoaxs = 101325*Ma/(Za*R*288.15)           
        	rhows  = 1333*Mw/(Zw*R*293.15)             
        	rhoa   = self.p(x)*Ma/(self.Z(x,T,xw)*R*T)*(1-xw)       
        	rhow   = self.p(x)*Mw/(self.Z(x,T,xw)*R*T)*xw         
        	nprop = 1 + ((rhoa/rhoaxs)*(naxs-1) + (rhow/rhows)*(nws-1))
        	return nprop

#writing homogenous deferential equation in terms of y(gamma=lorentz factor) and x(height)----------------------------------
#dydx=dy/dx=d(gamma)/dx(here y=gamma) 
	def dydx(self, y, x0):
        	if y<1.:
            		result = 0.
        	else:
            		gammasq= np.square(y)
            		betasq = (1.-1./gammasq)

            		T1 = 2.*self.En*(betasq)*gammasq/(2.+2.*y)
            		rho = self.rho0*np.exp(-x0/self.h0)
			Isq = np.square(10.*self.at_no*1e-12)
            		f1 = (self.K*self.zsq*self.at_no*rho/(self.Am*betasq*self.En))
            		f2 = ((0.5*np.log(2.*self.En*betasq*gammasq*T1/Isq)) - betasq)
			f3 = np.sqrt(gammasq)/self.Rad_len
            		result = (f1*f2)+f3 
        	return result

#solving the homogenous ordinary differential equation in terms of y and x or gamma and x-----------------------------------------
# calculate here the 'Gamma(lorentz factor)' of the electron at height 'h' for various initial conditions of 'y(y0,y2,y3,y4,y5,y6,y7)'.
	def gamma_lorentz(self,x):
		self.sol1 = odeint(self.dydx, self.y1, x)
		self.sol2 = odeint(self.dydx, self.y2, x)
		self.sol3 = odeint(self.dydx, self.y3, x) 
		with open('/n1/retnakah/corsika-76400/run/gamma_Ioni.txt','w') as fo:
			for d in self.sol1:
				fo.write(str(d) + '\n')
		return self.sol1, self.sol2, self.sol3
	

#beta using gamma(beta=sqrt(1-1/gamma**2))-------------------------------------------------
	def beta_electron(self):
		self.sol1_beta1 = np.sqrt(1.-1./np.square(self.sol1))
		self.sol2_beta2 = np.sqrt(1.-1./np.square(self.sol2))
		self.sol3_beta3 = np.sqrt(1.-1./np.square(self.sol3))		
		self.sol1_beta1[np.isnan(self.sol1_beta1)] = 1e-100
		with open('/n1/retnakah/corsika-76400/run/trial_beta.txt','w') as fo:
				for d in self.sol1_beta1:
					fo.write(str(d) + '\n')
		
		return self.sol1_beta1, self.sol2_beta2, self.sol3_beta3 

#calculating minimum gamma for minimum(300nm) and maximum wavelength(700nm) using refractive index calculated above:----------------
#eta1 = refractive index(n)-1
#gamma1=(1+eta)/(eta*(sqrt(1+2/eta)))
	def eta_gamma(self,x,l,l1):
		self.eta1 = self.n(l,x,15,0,450) - 1.#l=0.3
		self.eta2 = self.n(l1,x,15,0,450) - 1.#l1=0.7
		self.gamma1=(1+self.eta1)/(self.eta1*(np.sqrt(1+2/self.eta1)))
		self.gamma2 = (1.+self.eta2)/self.eta2/np.sqrt(1.+2./self.eta2) 
		beta1 = np.sqrt(1.-1./self.gamma1**2)#----------------------corresponding beta:
		beta2 = np.sqrt(1.-1./self.gamma2**2)
		return self.gamma1, self.gamma2


#method to calculate refractive index for blue and red light using corresponding eta values-------------------------------------------------------------
    	def nref_b_r(self):
        	self.nref_b = self.eta1 + 1 
        	self.nref_r = self.eta2 + 1
        	return self.nref_b, self.nref_r

#converting array of self.nref_b(for blue) values to a list, also converting array of multiple lists of self.sol1_beta1(\beta corresponding to \gamma0 = 120) in to a list-----------------
    	def conversion_blue(self):
        	tel1 = self.nref_b
        	self.tel2 = tel1.tolist()#-------------------converting array of list (refractive index corresponding to \lambda=300) to a list 
        	self.tel3 = []           #----------------------- converting array of multiple lists of self.sol1_beta1(\beta corresponding to \gamma0 = 120) in to a list 
        	for i in self.sol1_beta1:
            		self.tel3.append(i[0])
        	return self.tel3, self.tel2

	def conversion_1500(self):
		self.tel3_1500 = []
		for j in self.sol2_beta2:
			self.tel3_1500.append(j[0])
		return self.tel3_1500

	def conversion_10000(self):
		self.tel3_10000 = []
		for k in self.sol3_beta3:
			self.tel3_10000.append(k[0])
		return self.tel3_10000
		
#converting array of self.nref(for red) values to a list, also converting array of multiple lists of self.sol1_beta1 in to a list----------------------- 
    	def conversion_red(self):
        	tel11 = self.nref_r
        	self.tel4 = tel11.tolist() #-------------------------converting array of list (refractive index corresponding to \lambda=700) to a list 
        	return self.tel4

#calculating 'beta*refractive index'(both 300 and 700 nm) for all \gamma0 values then the "angle of emission"----------------------------------------
#only for one value of \gamma0=120 is given here: change for different \gamma0 values = 70, 80
    	def beta_n(self):
        	self.bn3 = np.multiply(self.tel3,self.tel2)#---------------------------for blue(\gamma0 = 100)
        	self.theta3 = 180.0/math.pi*np.arccos(1./self.bn3)
        	self.theta3[np.isnan(self.theta3)] = 1e-100
        	self.bn4 = np.multiply(self.tel4,self.tel3) #---------------------------for red(\gamma0 = 100)
		self.theta4 = 180.0/math.pi*np.arccos(1./self.bn4)
		self.theta4[np.isnan(self.theta4)] = 1e-100
		self.bn3_1500 = np.multiply(self.tel3_1500,self.tel2)#---------------------------for blue(\gamma0 = 1500)
        	self.theta3_1500 = 180.0/math.pi*np.arccos(1./self.bn3_1500)
		self.bn4_1500 = np.multiply(self.tel4,self.tel3_1500) #---------------------------for red(\gamma0 = 1500)
		self.theta4_1500 = 180.0/math.pi*np.arccos(1./self.bn4_1500)
		self.bn3_10000 = np.multiply(self.tel3_10000,self.tel2)#---------------------------for blue(\gamma0 = 10000)
        	self.theta3_10000 = 180.0/math.pi*np.arccos(1./self.bn3_10000)
		self.bn4_10000 = np.multiply(self.tel4,self.tel3_10000) #---------------------------for red(\gamma0 = 10000)
		self.theta4_10000 = 180.0/math.pi*np.arccos(1./self.bn4_10000)
        	return self.theta3_10000, self.theta4_1500, self.theta3_1500, self.theta4_1500, self.theta3, self.theta4

#calculating the radius of the electron rings obtained using radius = tan(theta)*height--------------------------------------------------------
#only for one value of \gamma0=120 is given here: change for different \gamma0 values = 70, 80
	def radius(self,x):
		theta_rad_b = (math.pi/180.0)*self.theta3#----------------------------converting angle corresponding to ring of blue light from degree to radians(100)
		theta_rad_r = (math.pi/180.0)*self.theta4 #---------------------------converting angle corresponding to ring of red light from degree to radians(100)
		theta_rad_b_1500 = (math.pi/180.0)*self.theta3_1500#----------------------------converting angle corresponding to ring of blue light from degree to radians(1500)
		theta_rad_r_1500 = (math.pi/180.0)*self.theta4_1500 #---------------------------converting angle corresponding to ring of red light from degree to radians(1500)
		theta_rad_b_10000 = (math.pi/180.0)*self.theta3_10000#----------------------------converting angle corresponding to ring of blue light from degree to radians(10000)
		theta_rad_r_10000 = (math.pi/180.0)*self.theta4_10000#---------------------------converting angle corresponding to ring of red light from degree to radians(10000)
		self.radius_b = x*np.tan(theta_rad_b)#-----------------------------------------for blue(100)
		self.radius_r = x*np.tan(theta_rad_r)#----------------------------------------for red(100)
		self.radius_b_1500 = x*np.tan(theta_rad_b_1500)#-----------------------------------------for blue(1500)
		self.radius_r_1500 = x*np.tan(theta_rad_r_1500)#----------------------------------------for red(1500)
		self.radius_b_10000 = x*np.tan(theta_rad_b_10000)#-----------------------------------------for blue(10000)
		self.radius_r_10000 = x*np.tan(theta_rad_r_10000)#----------------------------------------for red(10000)
		return self.radius_b, self.radius_r, self.radius_b_1500, self.radius_r_1500, self.radius_b_10000, self.radius_r_10000

#calculating the number of photons emitted per distance travelled(i.e. dN/dX with x )-------------------------------------------------------------
#Equation for dN/dx is taken from Particle Data Group(PDG)------------------------100 
	def dndx_num(self, x):
        	self.fbeta = interpolate.interp1d(x, self.tel3, kind='linear')
		self.Result3 = []		
		for j in self.fbeta(x):
			func = lambda l:((1-1/((j**2)*(self.n(l,x,15,0,450)[j])**2)))*(2*np.pi*self.alpha0*self.zsq*10**6/l**2)#------------------10**6 is added to convert one of the \lambda which would be remaining from lambda**2 after the integration in micrometer to meter to get the end value in "number per meter"(other lambda will be cancelled by the lambda in the integration limit)
        		result3, _ = quad(func, 0.3, 0.7)
			self.Result3.append(result3)	
		return self.Result3

#--------------------------------------1500
	def dndx_num_1500(self, x):
        	self.fbeta_1500 = interpolate.interp1d(x, self.tel3_1500, kind='linear')
		self.Result3_1500 = []		
		for k in self.fbeta_1500(x):
			func_1500 = lambda l:((1-1/((k**2)*(self.n(l,x,15,0,450)[k])**2)))*(2*np.pi*self.alpha0*self.zsq*10**6/l**2)#------------------10**6 is added to convert one of the \lambda which would be remaining from lambda**2 after the integration in micrometer to meter to get the end value in "number per meter"(other lambda will be cancelled by the lambda in the integration limit)
        		result3_1500, _ = quad(func_1500, 0.3, 0.7)
			self.Result3_1500.append(result3_1500)	
		return self.Result3_1500
#------------------------------------------------10000
	def dndx_num_10000(self, x):
        	self.fbeta_10000 = interpolate.interp1d(x, self.tel3_10000, kind='linear')
		self.Result3_10000 = []		
		for i in self.fbeta_10000(x):
			func_10000 = lambda l:((1-1/((i**2)*(self.n(l,x,15,0,450)[i])**2)))*(2*np.pi*self.alpha0*self.zsq*10**6/l**2)#------------------10**6 is added to convert one of the \lambda which would be remaining from lambda**2 after the integration in micrometer to meter to get the end value in "number per meter"(other lambda will be cancelled by the lambda in the integration limit)
        		result3_10000, _ = quad(func_10000, 0.3, 0.7)
			self.Result3_10000.append(result3_10000)	
		return self.Result3

#-------------------------------------------------------
#calculating the total number of photons(N) from dN/dx by integrating over x-----------------------------------------------N = 19992359.2906
	def N_dndx(self,x):
		N_dndx = trapz(self.Result3, (405.581116223,500))#----------------limits are made in to meters		
		return N_dndx

#TRIAL
#calculating angle of emission for wavelength(0.5 micrometer)--------------------------
	def theta(self,x):
		self.fbeta = interpolate.interp1d(x, self.tel3, kind='linear')
		#l = np.linspace(0.3,0.7,5000)
		self.theta0 = (np.arccos(1/(np.multiply(self.fbeta(x),self.n(0.5,x,15,0,450)))))#in radian
		self.theta0[np.isnan(self.theta0)] = 0.0
		with open('/n1/retnakah/corsika-76400/run/trial.txt','w') as fo:
				for d in self.theta0:
					fo.write(str(d) + '\n')
		return self.theta0

#calculating the slope of theta_ch vs x curve, d\theta/dx-----------------------------------------------
	def slope_dtheta_dx(self,x):
		self.func1 = (diff(self.theta0)/diff(x))*10**2 #-----------------#10**2 is to convert cm inverse to meter inverse
		return self.func1

#calculating the number of photons dN/d\Omega : Function is modified using the property of delta function(derivative)-------------------------------
#fixing \lambda=0.5 micro meter and e+6 is added to cancel one unit of length in \lambda in denomenator (convert micrometer to meter--->)end result will be in number per steradian per micro meter 
#0.1micro meter is multiplied to make the dimension correct to compensate to the \lambda in micrometer in the left side denomenator hence the output is in number per steradian
	def dndomega(self,x):
		self.Result4 = []
		for k in self.theta0:
			func2 = 0.000000006*np.sin(k)*(1/np.absolute(self.func1[k]))*(self.alpha0*0.1*self.zsq*10**6/0.5**2) 
			self.Result4.append(func2)
		with open('/n1/retnakah/corsika-76400/run/weights.txt','w') as fo:
				for d in self.Result4:
					fo.write(str(d) + '\n')
		return self.Result4
#no.per pixel
#	def no_pixel(self,x):
#		self.pixel = []
#		for l in self.Result4:
#			number_pixel = 2*np.pi(1-np.cos(self.theta_pixel_hess))*self.Result4
#			self.pixel.append(number_pixel)
#		return pixel

#comparing the results from analytical calculation to the results from CORSIKA simulation
#Converting the dN/d\Omega vs \theta_ch to a 2D histogram(u-v) to make direct comparison with the 2D histogram(u-v) from CORSIKA simulation
#First we are making a grid of u and v and fill it with photons that we got from dN/d\Omega calculation and finally we have to normalise it using the solid angle to do a solid comparison between the two.
#u=sin(\theta)*cos(\phi) and v=sin(\theta)*sin(\phi). Here \theta=\theta_ch(obtained from the above calculation) 
	def u_v_conversion(self,phi1):
		self.u = np.sin(self.theta0)*np.cos(phi1)#2*k+l
		self.v = np.sin(self.theta0)*np.sin(phi1)#2*k+k**2 
		self.value = (180/np.pi)*np.arcsin(np.sqrt(self.u**2 + self.v**2)) #--------------calculating \theta_ch(here it is defined as 'value') back from u and v obtained abouve: theta_ch = arcsin(sqrt(u^2 + v^2))
		with open('/n1/retnakah/corsika-76400/run/v_brem.txt','w') as fo:
				for d in self.v:
					fo.write(str(d) + '\n')
		return self.u, self.v, self.value
'''
def u_v_conversion(self,l):
		self.l1 = l*math.pi/180	
		self.u1 = []
		self.v1 = []
		for k in self.theta0:
			u = np.sin(k)*np.cos(self.l1)#2*k+l
			v = np.sin(k)*np.sin(self.l1)#2*k+k**2 
			self.u1.append(u)
			self.v1.append(v)
		return self.u1, self.v1
'''




'''
#ORIGINAL
#calculating angle of emission for a range of wavelength(0.3-0.7 micrometer)--------------------------
	def theta(self,x):
		self.fbeta = interpolate.interp1d(x, self.tel3, kind='linear')
		l = np.linspace(0.3,0.7,5000)
		self.theta0 = (np.arccos(1/(np.multiply(self.fbeta(x),self.n(l,x,15,0,450)))))#in radian
		self.theta0[np.isnan(self.theta0)] = 1e-100
		with open('/n1/retnakah/corsika-76400/run/trial.txt','w') as fo:
				for d in self.theta0:
					fo.write(str(d) + '\n')
		return self.theta0
'''
#ORIGINAL
'''
#calculating the number of photons dN/d\Omega : Function is modified using the property of delta function(derivative)-------------------------------
#fixing \lambda=0.5 micro meter and e+6 is added to cancel one unit of length in \lambda in denomenator (convert micrometer to meter--->)end result will be in number per steradian per micro meter 
#0.1micro meter is multiplied to make the dimension correct to compensate to the \lambda in micrometer in the left side denomenator hence the output is in number per steradian
	def dndomega(self,x):
		self.Result4 = []
		for k in self.tel3:
			func2 = np.sin(np.arccos(1/(k*self.n(0.5,x*10**-2,15,0,450)[k])))*(1/np.absolute(self.func1[k]))*(self.alpha0*0.1*self.zsq*10**6/0.5**2) 
			self.Result4.append(func2)
		return self.Result4



def brem_diff(self, e, x):#---------------------------------------------------------------------------------------------------Bremstrahlung loss
		if e<1.:
			brem_result = 0.
		else:
			brem_result = e/self.Rad_len
		return brem_result

	def bremsstrahlung(self,x):
		self.brem1 = odeint(self.brem_diff, self.y1, x)
		self.brem2 = odeint(self.brem_diff, self.y2, x)
		self.brem3 = odeint(self.brem_diff, self.y3, x)
		#self.brem1 = self.y1*np.exp(-1*x/self.Rad_len)
		#self.brem2 = self.y2*np.exp(-1*x/self.Rad_len)
		#self.brem3 = self.y3*np.exp(-1*x/self.Rad_len)
		#with open('/n1/retnakah/corsika-76400/run/gamma_brem.txt','w') as fo:
		#	for d in self.brem1:
		#		fo.write(str(d) + '\n')
		return self.brem1, self.brem2, self.brem3

	def total_momentum_100(self,x):
		self.gamma1_new = []
		
		for i in range(0, len(self.sol1)):
			self.gamma1_new.append(self.sol1[i]+self.brem1[i])
		with open('/n1/retnakah/corsika-76400/run/gamma_new.txt','w') as fo:
			for d in self.gamma1_new:
				fo.write(str(d) + '\n')
		return self.gamma1_new

	def total_momentum_500(self,x):
		self.gamma2_new = []
		for j in range(0, len(self.sol2)):
			self.gamma2_new.append(self.sol2[j]+self.brem2[j])
		return self.gamma2_new

	def total_momentum_1000(self,x):
		self.gamma3_new = []
		for j in range(0, len(self.sol3)):
			self.gamma3_new.append(self.sol3[j]+self.brem3[j])
		return self.gamma3_new
'''
