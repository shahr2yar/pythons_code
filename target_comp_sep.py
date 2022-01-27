from sys import argv
from matplotlib.pyplot import *
from matplotlib.legend_handler import HandlerLine2D
from numpy import *
from math import *
import glob
from scipy import * 
import os

script,argv1,argv2=argv
mnt=int(argv1)	#mass or number of target star
ns=int(argv2)	#number of total barionic stars
ts=0.148786818949782E+20/(365.242199*24*3600*1e6) #time scale in megayear to G=1

#/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\
#/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\	
rtarget=[]
vtarget=[]
target_p=[]
target_mom=[]
target_angular=[]
target_angula=0

rtarget_cm=[]
vtarget_cm=[]
target_p_cm=[]
target_mom_cm=[]
target_angular_cm=[]


m_tot=0
x_cm=0
y_cm=0
z_cm=0
vx_cm=0
vy_cm=0
vz_cm=0
xcm=[]
ycm=[]
zcm=[]
vxcm=[]
vycm=[]
vzcm=[]

files=glob.glob("./md/data/list*.dat")
files=sorted(files)
nh=0
nt=ns+nh

for f in files:
	print(f)
	data=genfromtxt(f)
	close
	
	for i in range(nt):	#find total barionic mass
		m_tot=m_tot+data[i,6]
			
	for i in range(nt):	#find center of mass
		x_cm=x_cm+data[i,6]*data[i,0]/m_tot
		y_cm=y_cm+data[i,6]*data[i,1]/m_tot
		z_cm=z_cm+data[i,6]*data[i,2]/m_tot
		vx_cm=vx_cm+data[i,6]*data[i,3]/m_tot
		vy_cm=vy_cm+data[i,6]*data[i,4]/m_tot
		vz_cm=vx_cm+data[i,6]*data[i,5]/m_tot
	xcm.append(x_cm)
	ycm.append(y_cm)	
	zcm.append(z_cm)
	vxcm.append(vx_cm)
	vycm.append(vy_cm)
	vzcm.append(vz_cm)
	
	for i in range(nt):	#determine heavy star position, velocity and angular momentum
		if data[i,7]==mnt:
			#target star distance 
			rtarget.append(sqrt(data[i,0]**2+data[i,1]**2+data[i,2]**2))	
			rtarget_cm.append(sqrt((data[i,0]-xcm[-1])**2+(data[i,1]-ycm[-1])**2+(data[i,2]-zcm[-1])**2))
			#magnitude of target star velocity 
			vtarget.append(sqrt(data[i,3]**2+data[i,4]**2+data[i,5]**2))	
			vtarget_cm.append(sqrt((data[i,3]-vxcm[-1])**2+(data[i,4]-vycm[-1])**2+(data[i,5]-vzcm[-1])**2))	
			#target star position vector
			target_p.append([data[i,0],data[i,1],data[i,2]])		
			target_p_cm.append([data[i,0]-xcm[-1],data[i,1]-ycm[-1],data[i,2]-zcm[-1]])
			#momentum of target star		
			mt=data[i,6]			
			target_mom=[(data[i,3])*mt,(data[i,4])*mt,(data[i,5])*mt]				
			target_mom_cm=[(data[i,3]-vxcm[-1])*mt,(data[i,4]-vycm[-1])*mt,(data[i,5]-vzcm[-1])*mt]	
			#angular momentum of target star 
			import numpy as np
			target_angula=cross(target_p[-1],target_mom)
			target_angular.append(np.linalg.norm(target_angula))
			target_angula_cm=cross(target_p_cm[-1],target_mom_cm)
			target_angular_cm.append(np.linalg.norm(target_angula_cm))	
	m_tot=0
	x_cm=0
	y_cm=0
	z_cm=0
	vx_cm=0
	vy_cm=0
	vz_cm=0
		
tar_pos11=[]
tar_ang12=[]
for i in range(len(files)):
	tar_pos11.append(rtarget_cm[i]/rtarget_cm[0])
	tar_ang12.append(target_angular_cm[i]/target_angular_cm[0])
for i in range(len(tar_pos11)):
	if(abs(tar_pos11[i]-0.25)<0.01):
		im1=i
		break
#/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\
rtarget=[]
vtarget=[]
target_p=[]
target_mom=[]
target_angular=[]
target_angula=0

rtarget_cm=[]
vtarget_cm=[]
target_p_cm=[]
target_mom_cm=[]
target_angular_cm=[]


m_tot=0
x_cm=0
y_cm=0
z_cm=0
vx_cm=0
vy_cm=0
vz_cm=0
xcm=[]
ycm=[]
zcm=[]
vxcm=[]
vycm=[]
vzcm=[]

files=glob.glob("./nd/data/list*.dat")
files=sorted(files)
nh=100000
nt=ns+nh

for f in files:
	print(f)
	data=genfromtxt(f)
	close
	
	for i in range(nt):	#find total barionic mass
		if data[i,6]==8.40892E+05:
			continue
		else:
			m_tot=m_tot+data[i,6]

			
	for i in range(nt):	#find center of mass
		if data[i,6]==8.40892E+05:
			continue
		else:
			x_cm=x_cm+data[i,6]*data[i,0]/m_tot
			y_cm=y_cm+data[i,6]*data[i,1]/m_tot
			z_cm=z_cm+data[i,6]*data[i,2]/m_tot
			vx_cm=vx_cm+data[i,6]*data[i,3]/m_tot
			vy_cm=vy_cm+data[i,6]*data[i,4]/m_tot
			vz_cm=vx_cm+data[i,6]*data[i,5]/m_tot
	xcm.append(x_cm)
	ycm.append(y_cm)	
	zcm.append(z_cm)
	vxcm.append(vx_cm)
	vycm.append(vy_cm)
	vzcm.append(vz_cm)
	for i in range(nt):	#determine heavy star position, velocity and angular momentum
		if data[i,7]==mnt:
			#target star distance 
			rtarget.append(sqrt(data[i,0]**2+data[i,1]**2+data[i,2]**2))	
			rtarget_cm.append(sqrt((data[i,0]-xcm[-1])**2+(data[i,1]-ycm[-1])**2+(data[i,2]-zcm[-1])**2))
			#magnitude of target star velocity 
			vtarget.append(sqrt(data[i,3]**2+data[i,4]**2+data[i,5]**2))	
			vtarget_cm.append(sqrt((data[i,3]-vxcm[-1])**2+(data[i,4]-vycm[-1])**2+(data[i,5]-vzcm[-1])**2))	
			#target star position vector
			target_p.append([data[i,0],data[i,1],data[i,2]])		
			target_p_cm.append([data[i,0]-xcm[-1],data[i,1]-ycm[-1],data[i,2]-zcm[-1]])
			#momentum of target star		
			mt=data[i,6]			
			target_mom=[(data[i,3])*mt,(data[i,4])*mt,(data[i,5])*mt]				
			target_mom_cm=[(data[i,3]-vxcm[-1])*mt,(data[i,4]-vycm[-1])*mt,(data[i,5]-vzcm[-1])*mt]	
			#angular momentum of target star 
			import numpy as np
			target_angula=cross(target_p[-1],target_mom)
			target_angular.append(np.linalg.norm(target_angula))
			target_angula_cm=cross(target_p_cm[-1],target_mom_cm)
			target_angular_cm.append(np.linalg.norm(target_angula_cm))	
	print(n_test)
	m_tot=0
	x_cm=0
	y_cm=0
	z_cm=0
	vx_cm=0
	vy_cm=0
	vz_cm=0
		
tar_pos21=[]
tar_ang22=[]
for i in range(len(files)):
	tar_pos21.append(rtarget_cm[i]/rtarget_cm[0])
	tar_ang22.append(target_angular_cm[i]/target_angular_cm[0])
for i in range(len(tar_pos21)):
	if(abs(tar_pos21[i]-0.25)<0.01):
		im2=i
		break	
#/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\
rtarget=[]
vtarget=[]
target_p=[]
target_mom=[]
target_angular=[]
target_angula=0

rtarget_cm=[]
vtarget_cm=[]
target_p_cm=[]
target_mom_cm=[]
target_angular_cm=[]


m_tot=0
x_cm=0
y_cm=0
z_cm=0
vx_cm=0
vy_cm=0
vz_cm=0
xcm=[]
ycm=[]
zcm=[]
vxcm=[]
vycm=[]
vzcm=[]

files=glob.glob("./pure_nd/data/list*.dat")
files=sorted(files)
nh=0
nt=ns+nh
for f in files:
	print(f)
	data=genfromtxt(f)
	close
	
	for i in range(nt):	#find total barionic mass
		m_tot=m_tot+data[i,6]
			
	for i in range(nt):	#find center of mass
		x_cm=x_cm+data[i,6]*data[i,0]/m_tot
		y_cm=y_cm+data[i,6]*data[i,1]/m_tot
		z_cm=z_cm+data[i,6]*data[i,2]/m_tot
		vx_cm=vx_cm+data[i,6]*data[i,3]/m_tot
		vy_cm=vy_cm+data[i,6]*data[i,4]/m_tot
		vz_cm=vx_cm+data[i,6]*data[i,5]/m_tot
	xcm.append(x_cm)
	ycm.append(y_cm)	
	zcm.append(z_cm)
	vxcm.append(vx_cm)
	vycm.append(vy_cm)
	vzcm.append(vz_cm)
	for i in range(nt):	#determine heavy star position, velocity and angular momentum
		if data[i,7]==mnt:
			#target star distance 
			rtarget.append(sqrt(data[i,0]**2+data[i,1]**2+data[i,2]**2))	
			rtarget_cm.append(sqrt((data[i,0]-xcm[-1])**2+(data[i,1]-ycm[-1])**2+(data[i,2]-zcm[-1])**2))
			#magnitude of target star velocity 
			vtarget.append(sqrt(data[i,3]**2+data[i,4]**2+data[i,5]**2))	
			vtarget_cm.append(sqrt((data[i,3]-vxcm[-1])**2+(data[i,4]-vycm[-1])**2+(data[i,5]-vzcm[-1])**2))	
			#target star position vector
			target_p.append([data[i,0],data[i,1],data[i,2]])		
			target_p_cm.append([data[i,0]-xcm[-1],data[i,1]-ycm[-1],data[i,2]-zcm[-1]])
			#momentum of target star		
			mt=data[i,6]			
			target_mom=[(data[i,3])*mt,(data[i,4])*mt,(data[i,5])*mt]				
			target_mom_cm=[(data[i,3]-vxcm[-1])*mt,(data[i,4]-vycm[-1])*mt,(data[i,5]-vzcm[-1])*mt]	
			#angular momentum of target star 
			import numpy as np
			target_angula=cross(target_p[-1],target_mom)
			target_angular.append(np.linalg.norm(target_angula))
			target_angula_cm=cross(target_p_cm[-1],target_mom_cm)
			target_angular_cm.append(np.linalg.norm(target_angula_cm))	
	m_tot=0
	x_cm=0
	y_cm=0
	z_cm=0
	vx_cm=0
	vy_cm=0
	vz_cm=0
		
tar_pos31=[]
tar_ang32=[]
for i in range(len(files)):
	tar_pos31.append(rtarget_cm[i]/rtarget_cm[0])
	tar_ang32.append(target_angular_cm[i]/target_angular_cm[0])
for i in range(len(tar_pos31)):
	if(abs(tar_pos31[i]-0.25)<0.01):
		im3=i
		break
#/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\
#/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\

t1=[]
f=open("./md/data/time_step","r");
for line in f:
	t1.append(float(line)*ts)	
tm1=t1[im1]*ts
t2=[]
f=open("./nd/data/time_step","r");
for line in f:
	t2.append(float(line)*ts)	
tm2=t2[im2]*ts
t3=[]
f=open("./pure_nd/data/time_step","r");
for line in f:
	t3.append(float(line)*ts)	
tm3=t3[im3]*ts

############################################################################
############################################################################

figure(figsize=(14,7),dpi=500)

subplot(311)
plot(t1[0:len(t1):1],tar_ang12[0:len(t1):1],linewidth=1,label='Milgromian')
plot(t3[0:len(t3):1],tar_ang32[0:len(t3):1],linewidth=1,linestyle='--',label="Newtonian\nwithout DM")
yticks([1.0,0.75,0.5,0.25,0.0],fontsize=14)
xticks(fontsize=15)
ylabel(r'$\frac{L}{L_0}$',fontsize=20)
legend(loc='upper right',fontsize=11)
text_pos_x =-30
text_pos_y =0.6
text(text_pos_x, text_pos_y, r"$r(t=0)=r_h=15\;kpc$",fontsize=12)
text_pos_x =-30
text_pos_y =0.4
text(text_pos_x, text_pos_y, r"$v(t=0)=v_c$",fontsize=14)
text_pos_x =-30
text_pos_y =0.2
text(text_pos_x, text_pos_y,r"$\frac{a_n(t=0)}{a_0}\sim 0.003$",fontsize=12)

subplot(312)
plot(t1[0:len(t1):1],tar_ang12[0:len(t1):1],linewidth=1,label='Milgromian')
plot(t2[0:len(t2):1],tar_ang22[0:len(t2):1],linewidth=1,label="Newtonian\nwith DM")
yticks([1.0,0.75,0.5,0.25,0.0],fontsize=15)
legend(loc='upper right',fontsize=11)
xticks(fontsize=15)
ylabel(r'$\frac{L}{L_0}$',fontsize=20)
xlabel(r'$t\;(Myr)$',fontsize=20)

savefig("target_angular_evolution_wrt_cm.eps")
savefig("target_angular_evolution_wrt_cm.jpg")
close()

#/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\
#/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\/|\

figure(figsize=(14,7),dpi=500)

subplot(211)
plot(t1[0:len(t1):1],tar_pos11[0:len(t1):1],linewidth=1,label='Milgromian')
plot(t3[0:len(t3):1],tar_pos31[0:len(t3):1],linewidth=1,linestyle='--',label="Newtonian\nwithout DM")
yticks([1.0,0.75,0.5,0.25,0.0],fontsize=14)
xticks(fontsize=15)
ylabel(r'$\frac{D}{D_0}$',fontsize=20)
legend(loc='upper right',fontsize=14)
text_pos_x =-30
text_pos_y =0.6
text(text_pos_x, text_pos_y, r"$r(t=0)=r_h=15\;KPC$",fontsize=12)
text_pos_x =-30
text_pos_y =0.4
text(text_pos_x, text_pos_y, r"$v(t=0)=v_c$",fontsize=14)
text_pos_x =-30
text_pos_y =0.2
text(text_pos_x, text_pos_y,r"$\frac{a_n(t=0)}{a_0}\sim 0.003$",fontsize=12)

subplot(212)
plot(t1[0:len(t1):1],tar_pos11[0:len(t1):1],linewidth=1,label='Milgromian')
plot(t2[0:len(t2):5],tar_pos21[0:len(t2):5],linewidth=1,linestyle='--',label="Newtonian\nwith DM")
yticks([1.0,0.75,0.5,0.25,0.0],fontsize=15)
legend(loc='upper right',fontsize=14)
xticks(fontsize=15)
ylabel(r'$\frac{D}{D_0}$',fontsize=20)
xlabel(r'$t\;(Myr)$',fontsize=20)

savefig("target_distance_evolution_wrt_cm.eps")
savefig("target_distance_evolution_wrt_cm.jpg")
close()

############################################################################
############################################################################
print('-'*50)
