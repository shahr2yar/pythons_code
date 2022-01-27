import os
from sys import *
time=[]
for i in range(1,181):
	if i<10:
		os.chdir('./output_0000%d'%i)
		f=open("info_0000%d.txt"%i)
		j=0
		k=0
		for line in f:
			j+=1
			if j==9:
				for word in line.split( ):
					k+=1
					if k==3:
						time.append(word)
		os.chdir('../')
	elif i<100:
		os.chdir('./output_000%d'%i)
		f=open("info_000%d.txt"%i)
		j=0
		k=0
		for line in f:
			j+=1
			if j==9:
				for word in line.split( ):
					k+=1
					if k==3:
						time.append(word)
		os.chdir('../')
	elif i<1000:
		os.chdir('./output_00%d'%i)
		f=open("info_00%d.txt"%i)
		j=0
		k=0
		for line in f:
			j+=1
			if j==9:
				for word in line.split( ):
					k+=1
					if k==3:
						time.append(word)
		os.chdir('../')
	elif i<10000:
		os.chdir('./output_0%d'%i)
		f=open("info_0%d.txt"%i)
		j=0
		k=0
		for line in f:
			j+=1
			if j==9:
				for word in line.split( ):
					k+=1
					if k==3:
						time.append(word)
		os.chdir('../')




g=open("time_step","w")
for t in time:
	g.write(t)
	g.write('\n')

g.close()
