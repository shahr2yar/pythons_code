from sys import argv
from matplotlib.pyplot import *
from matplotlib.legend_handler import HandlerLine2D
from numpy import *
from math import *
import glob 
import os

script,argv1,argv2=argv
mnt=int(argv1)	#mass or number of target star
ns=int(argv2)	#number of total stars

m_b=2.01578E+06 #nd
#m_b=2.00345E+06 #md
ts=0.148786818949782E+20/(365.242199*24*3600*1e6) #time scale in megayear to G=1
files=glob.glob("*.dat")
files=sorted(files)
star_px=[]
star_py=[]
star_pz=[]
dm_px=[]
dm_py=[]
dm_pz=[]
target_px=[]
target_py=[]
target_pz=[]
k=0
t=[]
f=open("time_step","r");
for line in f:
	t.append(float(line)*ts)
for f in files:
	print(f)
	data=genfromtxt(f)
	close()
	
	for i in range(ns):	#determine all stars position, velocity and angular momentum
		if data[i,6]==m_b:
			star_px.append(data[i,0])		#all stars position vector
			star_py.append(data[i,1])
			star_pz.append(data[i,2])
		else:
			dm_px.append(data[i,0])
			dm_py.append(data[i,1])
			dm_pz.append(data[i,2])
		print(i)
	for i in range(ns):
		if data[i,7]==mnt:
			target_px.append(data[i,0])		#target star position vector
			target_py.append(data[i,1])
			target_pz.append(data[i,2])
	subplot(211)
	scatter(dm_px,dm_py,c='blue',s=0.1)
	scatter(star_px,star_py,marker='*',s=1)
	scatter(target_px,target_py,c='red',marker='o',s=1)
	tt=float(t[k])
	title(r'$t=%0.2f\;Myr$'%tt,fontsize=14 )
	k+=1
	ylabel('Y (kpc)')
	axis('equal')
	subplot(212)
	scatter(dm_px,dm_pz,c='blue',s=0.1)
	scatter(star_px,star_pz,marker='*',s=1)
	scatter(target_py,target_pz,c='red',marker='o',s=1)
	xlabel('X (kpc)')
	ylabel('Z (kpc)')
	axis('equal')
	savefig("galaxy evolution in %s snapshot.png"%f)
	close()
	star_px=[]
	star_py=[]
	star_pz=[]
	dm_px=[]
	dm_py=[]
	dm_pz=[]
	print('*'*10)
############################################################################
print('-'*50)
