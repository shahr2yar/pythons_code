from numpy import *
from math import *
import glob
#from scipy import * 
from sys import argv

script, r1, r2, r3=argv #r=[r1,r2,r3] position of target star, rh is mid mass radius
r=[]
r.append(float(r1))
r.append(float(r2))
r.append(float(r3))

m_b=10251.474427
rr=sqrt(r[0]**2+r[1]**2+r[2]**2)
print('hi\n')
#physical constant
#g=6.67408e-11 #m**3/(kg*s**2)
#g=0.004498 #pc**3/(m_sun*Myr**2) or kpc**3/(1e9*m_sun*Myr**2)
#a0=1.12*10**(-10) #m/s**2
#a0=0.00361	#kpc/myr^2

masst_b=0;masst=0;massr=0;massr_b=0
nt_b=0
nr_b=0
nt=0
files=glob.glob("*.dat")
for f in files:
	data=genfromtxt(f)
	for i in range(len(data)):
		if int(data[i,6])==int(m_b):
			masst_b+=data[i,6]
			nt_b+=1
		masst+=data[i,6]
		nt+=1	
		if((data[i,0]**2+data[i,1]**2+data[i,2]**2)**0.5<=rr):
			massr+=data[i,6]
			if(int(data[i,6])==int(m_b)):
				massr_b+=data[i,6]
				nr_b+=1
 
m_dm=str(masst-masst_b)
m_ratio=float(m_dm)/masst_b
from math import *
mratio=100*massr/masst
print("target in selected radious is above %f"%mratio+' percent of total mass')
g=0.004492e-9	#kpc**3/(m_sun*myr**2)
a0=0.00361	#kpc/myr^2 
chdim=3.103e-8
a_t=g*massr/rr**2	#kpc/myr**2
if a_t<a0:
	print("your galaxy with acceleration in your selected radius, %f (kpc/myr^2) is in milgromian regime"%a_t)
else:
	print("your galaxy with acceleration in your selected radius, %f (kpc/myr^2) is in newtonian regime"%a_t)


from math import *
class vf(object):
	def __init__(self,r,mt,mb):
		self.r=r
		self.mt=mt
		self.mb=mb
	def newtonian(self):	#derive newtonian circular velocity
		g=0.004492e-9	#kpc**3/(m_sun*myr**2)
		a0=0.00361	#kpc/myr^2  
		v_n=sqrt(g*self.mt/self.r)	#kpc/Myr
		return v_n
	def milgromian(self):	#derive milgromian circular velocity
		#  nu = 0.5 + 0.5*sqrt(x*x + 4.0*x)/x              Simple nu function
		#  nu = sqrt( 0.5 + sqrt( 0.25 + 1.0/x**2 ) ) 	   Standard nu function
		#  nu = 1.0/(1.0 - exp(-sqrt(x))) 
		g=0.004492e-9	#kpc**3/(m_sun*myr**2)
		a0=0.00361	#kpc/myr^2	
		v_n=sqrt(g*self.mb/self.r)	
		a_n=v_n**2/self.r
		x=a_n/a0
		a_m=(0.5+0.5*sqrt(x*x+4.0*x)/x)*a_n
		v_m=sqrt(self.r*a_m)
		return v_m


velof=vf(rr,massr,massr_b)
print("\n\n\nnewtonian velocity for cicular motion in effect of above mass is:")
print(velof.newtonian()*0.9785*1000)	#km/s
print("\n\n\nmilgromian velocity for cicular motion in effect of above mass is:")
print(velof.milgromian()*0.9785*1000)
g=0.004492e-9	#kpc**3/(m_sun*myr**2)
a0=0.00361	#kpc/myr^2
ic=open('IC','w')

if a_t<a0:
	ic.write("galaxy with a(r)= %f (kpc/Myr^2) < 0.00361 (kpc/Myr^2) is in the milgromian regime"%a_t)
else:
	ic.write("galaxy with a(r)= %f (kpc/Myr^2) > 0.00361 (kpc/Myr^2) is in the newtonian regime"%a_t)
velon=velof.newtonian()*0.9785*1000	#km/s
velom=velof.milgromian()*0.9785*1000
masst=float(masst)/1e9
massr=float(massr)/1e9
m_dm=float(m_dm)/1e9
masst_b=float(masst_b)/1e9
ic.write("\n\n\ntotal mass in this simulation is=%0.3f"%masst+"*1e9 solar mass")
ic.write("\n\n\ntotal mass underfoot of target is=%0.3f"%massr+"*1e9 solar mass")
ic.write("\n\n\nmass of dark matter in this simulation is=%0.3f"%m_dm+"*1e9 solar mass")
ic.write("\n\n\nmass of barionic matter in this simulation is=%0.3f "%masst_b+'*1e9 solar mass')
ic.write("\n\n\nratio of dark matter to barionic matter is=%0.3f"%m_ratio)
ic.write("\n\n\ntarget in selected radious (r=[%s,%s,%s] kpc) is above %0.3f"%(r1,r2,r3,mratio)+' percent of total mass')
ic.write("\n\n\nnewtonian velocity (km/s) for cicular motion in effect of inferior mass is: %0.12f"%velon)
ic.write("\n\n\nmilgromian velocity (km/s) for cicular motion in effect of inferior mass is: %0.12f"%velom)
ic.close()

print("-"*50)
