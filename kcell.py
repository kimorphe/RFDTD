#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt



class IMG:
    def __init__(self):
        self.nshow=0
    def load(self,fname):
        fp=open(fname,"r");
        fp.readline(); # computational domain
        self.time=float(fp.readline()); # computational domain
        fp.readline(); # computational domain
        Xa=np.array(fp.readline().lstrip().split(","));
        Xa=Xa.astype("float")

        fp.readline(); # computational domain
        Xb=np.array(fp.readline().lstrip().split(","));
        Xb=Xb.astype("float")

        fp.readline(); # computational domain
        Ndiv=np.array(fp.readline().lstrip().split(","));
        Ndiv=Ndiv.astype("int")
        fp.readline();	# Imaging area
        
        K=[];
        for row in fp:
            dat=row.lstrip().split(",");
            dat=list(map(float,dat))
            #K.append(dat[2]);
            K.append(dat[0]*dat[0]+dat[2]*dat[2]);
        fp.close();

        K=np.sqrt(np.array(K))
        #K=np.transpose(np.reshape(K,[Ndiv[0],Ndiv[1]]))
        K=np.transpose(np.reshape(K,[Ndiv[0],Ndiv[2]]))

        self.Xa=Xa;
        self.Xb=Xb;
        self.Ndiv=Ndiv;
        self.K=K

    def show(self,ax):
        #im=ax.imshow(K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",vmin=-0.01,vmax=0.01,cmap="jet");
        Xa=self.Xa;
        Xb=self.Xb;
        #im=ax.imshow(self.K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="bilinear",cmap="jet",vmin=-0.005,vmax=0.005);
        im=ax.imshow(self.K,origin="lower",extent=[Xa[0],Xb[0],Xa[2],Xb[2]],interpolation="bilinear",cmap="jet",vmin=0.0,vmax=0.003);
        if self.nshow==0:
            plt.colorbar(im);
        ax.set_aspect(1.0)
        ax.set_title("t="+str(self.time));
        self.nshow+=1;

if __name__=="__main__":
    fig=plt.figure();
    ax=fig.add_subplot(111)

    vdat=IMG()
    nums=np.arange(1,31,1)
    for k in nums:
        fname="v"+str(k)+".out";
        vdat.load(fname)
        vdat.show(ax)
        fout=fname.replace(".out",".png")
        print(fout)
        fig.savefig(fout,bbox_inches="tight")

