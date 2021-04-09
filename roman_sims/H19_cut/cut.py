import numpy as np
import glob

files=glob.glob('*norm_old_old.dat')
for fi in files:
	w,t=np.loadtxt(fi,unpack=True)
	#with open(fi[:-4]+'_old.dat','w') as f:
	#	for i in range(len(w)):
	#		f.write('%.1f %.20f\n'%(w[i],t[i]))
	t[t<1e-03]=0
	with open(fi[:-12]+'.dat','w') as f:
		for i in range(len(w)):
			f.write('%.1f %.20f\n'%(w[i],t[i]))
