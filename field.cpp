#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cij.h"
#include "field.h"

void Field::mem_alloc(int ndiv[3]){


	for(int i=0;i<3;i++){
		ngv[i]=ndiv[i]+1;	// number of grids for stress components 
		ngs[i]=ndiv[i]+2;	// number of grids for velocity components
	}

	S1=Field::mem_alloc3d(ngs);
	S2=Field::mem_alloc3d(ngs);
	S3=Field::mem_alloc3d(ngs);
	S4=Field::mem_alloc3d(ngs);
	S5=Field::mem_alloc3d(ngs);
	S6=Field::mem_alloc3d(ngs);

	V1=Field::mem_alloc3d(ngv);
	V2=Field::mem_alloc3d(ngv);
	V3=Field::mem_alloc3d(ngv);
};

double ***Field::mem_alloc3d(int ndiv[3]){

	Ndiv[0]=ndiv[0];
	Ndiv[1]=ndiv[1];
	Ndiv[2]=ndiv[2];
	int ndat=Ndiv[0]*Ndiv[1]*Ndiv[2];

	double *p1=(double *)malloc(sizeof(double)*ndat);
	double **p2=(double **)malloc(sizeof(double*)*Ndiv[0]*Ndiv[1]);
	double ***p3=(double ***)malloc(sizeof(double**)*Ndiv[0]);

	int i,j,k;
	for(k=0;k<ndat;k++) p1[k]=0.0;
	for(j=0;j<Ndiv[0]*Ndiv[1];j++) p2[j]=p1+j*Ndiv[2];
	for(i=0;i<Ndiv[0];i++) p3[i]=p2+i*Ndiv[1];
	
	int isum=0;
	for(i=0; i<Ndiv[0]; i++){
	for(j=0; j<Ndiv[1]; j++){
	for(k=0; k<Ndiv[2]; k++){
		p3[i][j][k]=isum++;
	}
	}
	}
	return(p3);
};

void Field::print_Vt(){
	int i,j,k;
	for(i=0; i<Ndiv[0]; i++){
	for(j=0; j<Ndiv[1]; j++){
	for(k=0; k<Ndiv[2]; k++){
		printf("%lf ",V2[i][j][k]);
	}
	printf("\n");
	}
	printf("\n");
	}
};
void Field::v2s(){
	int i,j,k;
	int i1,j1,k1;

	double phi[4];
	double dv11,dv12,dv13;
	double dv21,dv22,dv23;
	double dv31,dv32,dv33;

	for(i=0;i<ngs[0];i++){
	for(j=0;j<ngs[1];j++){
	for(k=0;k<ngs[2];k++){
	i1=i+1;
	j1=j+1;
	k1=k+1;

	phi[0]=V1[i1][j1][k1]-V1[i][j][k];
	phi[1]=V1[i1][j1][k] -V1[i][j][k1];
	phi[2]=V1[i1][j][k1] -V1[i][j1][k];
	phi[3]=V1[i1][j][k]  -V1[i][j1][k1];
	dv11=0.25*(phi[0]+phi[1]+phi[2]+phi[3]);
	dv12=0.25*(phi[0]+phi[1]-phi[2]-phi[3]);
	dv13=0.25*(phi[0]-phi[1]+phi[2]-phi[3]);

	phi[0]=V2[i1][j1][k1]-V2[i][j][k];
	phi[1]=V2[i1][j1][k] -V2[i][j][k1];
	phi[2]=V2[i1][j][k1] -V2[i][j1][k];
	phi[3]=V2[i1][j][k]  -V2[i][j1][k1];
	dv21=0.25*(phi[0]+phi[1]+phi[2]+phi[3]);
	dv22=0.25*(phi[0]+phi[1]-phi[2]-phi[3]);
	dv23=0.25*(phi[0]-phi[1]+phi[2]-phi[3]);

	phi[0]=V3[i1][j1][k1]-V3[i][j][k];
	phi[1]=V3[i1][j1][k] -V3[i][j][k1];
	phi[2]=V3[i1][j][k1] -V3[i][j1][k];
	phi[3]=V3[i1][j][k]  -V3[i][j1][k1];
	dv31=0.25*(phi[0]+phi[1]+phi[2]+phi[3]);
	dv32=0.25*(phi[0]+phi[1]-phi[2]-phi[3]);
	dv33=0.25*(phi[0]-phi[1]+phi[2]-phi[3]);
	}
	}
	}
};


void Field::s2v(){
	int i,j,k;
	int i1,j1,k1;

	double phi[4];
	double ds11,ds22,ds33;
	double ds42,ds43,ds51,ds53,ds61,ds62;

	for(i=0;i<ngv[0];i++){
	for(j=0;j<ngv[1];j++){
	for(k=0;k<ngv[2];k++){
	i1=i+1;
	j1=j+1;
	k1=k+1;
	phi[0]=S1[i1][j1][k1]-S1[i][j][k];
	phi[1]=S1[i1][j1][k] -S1[i][j][k1];
	phi[2]=S1[i1][j][k1] -S1[i][j1][k];
	phi[3]=S1[i1][j][k]  -S1[i][j1][k1];
	ds11=0.25*(phi[0]+phi[1]+phi[2]+phi[3]);

	phi[0]=S2[i1][j1][k1]-S2[i][j][k];
	phi[1]=S2[i1][j1][k] -S2[i][j][k1];
	phi[2]=S2[i1][j][k1] -S2[i][j1][k];
	phi[3]=S2[i1][j][k]  -S2[i][j1][k1];
	ds22=0.25*(phi[0]+phi[1]-phi[2]-phi[3]);

	phi[0]=S3[i1][j1][k1]-S3[i][j][k];
	phi[1]=S3[i1][j1][k] -S3[i][j][k1];
	phi[2]=S3[i1][j][k1] -S3[i][j1][k];
	phi[3]=S3[i1][j][k]  -S3[i][j1][k1];
	ds33=0.25*(phi[0]-phi[1]+phi[2]-phi[3]);

	phi[0]=S4[i1][j1][k1]-S4[i][j][k];
	phi[1]=S4[i1][j1][k] -S4[i][j][k1];
	phi[2]=S4[i1][j][k1] -S4[i][j1][k];
	phi[3]=S4[i1][j][k]  -S4[i][j1][k1];
	ds42=0.25*(phi[0]+phi[1]-phi[2]-phi[3]);
	ds43=0.25*(phi[0]-phi[1]+phi[2]-phi[3]);

	phi[0]=S5[i1][j1][k1]-S5[i][j][k];
	phi[1]=S5[i1][j1][k] -S5[i][j][k1];
	phi[2]=S5[i1][j][k1] -S5[i][j1][k];
	phi[3]=S5[i1][j][k]  -S5[i][j1][k1];
	ds51=0.25*(phi[0]+phi[1]+phi[2]+phi[3]);
	ds53=0.25*(phi[0]-phi[1]+phi[2]-phi[3]);

	phi[0]=S6[i1][j1][k1]-S6[i][j][k];
	phi[1]=S6[i1][j1][k] -S6[i][j][k1];
	phi[2]=S6[i1][j][k1] -S6[i][j1][k];
	phi[3]=S6[i1][j][k]  -S6[i][j1][k1];
	ds61=0.25*(phi[0]+phi[1]+phi[2]+phi[3]);
	ds62=0.25*(phi[0]+phi[1]-phi[2]-phi[3]);
	}
	}
	}
};

void DOMAIN::setup(
	double xa[3],
	double xb[3],
	double h
){
	dh=h;
	for(int i=0;i<3;i++){
		Xa[i]=xa[i];
		Xb[i]=xb[i];
		Ndiv[i]=ceil((Xb[i]-Xa[i])/dh);
		dx[i]=dh;
		Xb[i]=Xa[i]+dh*Ndiv[i];
	};

	printf("Xa=%lf %lf %lf\n",Xa[0],Xa[1],Xa[2]);
	printf("Xb=%lf %lf %lf\n",Xb[0],Xb[1],Xb[2]);
	printf("Ndiv=%d %d %d\n",Ndiv[0],Ndiv[1],Ndiv[2]);
	
};

