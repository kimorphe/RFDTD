#include<stdio.h>
#include<math.h>
#include<stdlib.h>


int main(){
	double PI=4.0*atan(1.0);
	int Ndiv[2]={100,100};
	int Ng[2];
	int i,j;

	double Xa[2]={ 0.0, 0.0};
	double Xb[2]={10.0,10.0};
	double Wd[2];
	double xx,yy;
	double dx[2];

	for(i=0;i<2;i++){
		Ng[i]=Ndiv[i]+1;
		Wd[i]=Xb[i]-Xa[i];
		dx[i]=Wd[i]/Ndiv[i];
	}


	FILE *fp=fopen("smp0.vtk","w");

	fprintf(fp,"# vtk DataFile Version 3.0\n");
	fprintf(fp,"test\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET STRUCTURED_GRID\n");
	fprintf(fp,"DIMENSIONS %d %d %d\n",Ng[0],Ng[1],1);
	fprintf(fp,"POINTS %d float\n",Ng[0]*Ng[1]);

	for(i=0;i<Ng[0];i++){
		xx=Xa[0]+dx[0]*i;
	for(j=0;j<Ng[1];j++){
		yy=Xa[1]+dx[1]*j;
		fprintf(fp," %lf %lf %lf\n",xx,yy,0.0);	
	}
	}

	double val;
	double kx,ky;

	fprintf(fp,"POINT_DATA %d\n",Ng[0]*Ng[1]);
	fprintf(fp,"SCALARS temp float\n");
	fprintf(fp,"LOOKUP_TABLE default\n");
	kx=2.*PI/(Wd[0]*0.2);
	ky=2.*PI/(Wd[1]*0.4);
	for(i=0;i<Ng[0];i++){
		xx=Xa[0]+dx[0]*i;
	for(j=0;j<Ng[1];j++){
		yy=Xa[1]+dx[1]*j;
		val=sin(kx*xx)*cos(ky*yy);
		fprintf(fp," %lf\n",val);	
	}
	}

	fclose(fp);


	return(0);
};
