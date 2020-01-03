#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt


fname="v40.out";

fp=open(fname,"r");

fp.readline(); # computational domain
Xa=np.array(fp.readline().lstrip().split(","));
Xa=Xa.astype("float")

fp.readline(); # computational domain
Xb=np.array(fp.readline().lstrip().split(","));
Xb=Xb.astype("float")

fp.readline(); # computational domain
Ndiv=np.array(fp.readline().lstrip().split(","));
Ndiv=Ndiv.astype("int")

print("Xa=",Xa)
print("Xb=",Xb)
print("Ndiv=",Ndiv)

fp.readline();	# Imaging area

K=[];
for row in fp:
    dat=row.lstrip().split(",");
    dat=list(map(float,dat))
    K.append(dat[2]);
    #K=np.hstack([K,np.array(dat)])
print(np.shape(K));
K=np.array(K)
K=np.transpose(np.reshape(K,[Ndiv[0],Ndiv[1]]))

fig=plt.figure();
ax=fig.add_subplot(111)
#im=ax.imshow(K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",vmin=-0.01,vmax=0.01,cmap="jet");
im=ax.imshow(K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",cmap="jet");
plt.colorbar(im);
plt.show();


fp.close();
