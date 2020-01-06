#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cij.h"
#include "field.h"

void Field::mem_alloc(int ndiv[3]){

	for(int i=0;i<3;i++){
		ngv[i]=ndiv[i];   // number of grids for velocity components 
		ngs[i]=ndiv[i]+1; // number of grids for stress   components
	}
//		Stress Grids
	S1=Field::mem_alloc3d(ngs);	// s11
	S2=Field::mem_alloc3d(ngs);	// s22
	S3=Field::mem_alloc3d(ngs);	// s33
	S4=Field::mem_alloc3d(ngs);	// s23
	S5=Field::mem_alloc3d(ngs);	// s31
	S6=Field::mem_alloc3d(ngs);	// s12
//		Velocity Grids
	V1=Field::mem_alloc3d(ngv);	// v1
	V2=Field::mem_alloc3d(ngv);	// v2
	V3=Field::mem_alloc3d(ngv);	// v3
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
		p3[i][j][k]=0.0;
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
	int i,j,k,I;
	int i1,j1,k1;

	double phi[4];
	double ds[6];
	double dv11,dv12,dv13;
	double dv21,dv22,dv23;
	double dv31,dv32,dv33;

	double dl=4.*dh/dt;


	for(i=1;i<ngs[0]-1;i++){
	for(j=1;j<ngs[1]-1;j++){
	for(k=1;k<ngs[2]-1;k++){
		i1=i-1;
		j1=j-1;
		k1=k-1;

		phi[0]=V1[i][j][k]   -V1[i1][j1][k1];
		phi[1]=V1[i][j][k1]  -V1[i1][j1][k];
		phi[2]=V1[i][j1][k]  -V1[i1][j][k1];
		phi[3]=V1[i][j1][k1] -V1[i1][j][k];
		dv11=(phi[0]+phi[1]+phi[2]+phi[3])/dl;
		dv12=(phi[0]+phi[1]-phi[2]-phi[3])/dl;
		dv13=(phi[0]-phi[1]+phi[2]-phi[3])/dl;

		phi[0]=V2[i][j][k]   -V2[i1][j1][k1];
		phi[1]=V2[i][j][k1]  -V2[i1][j1][k];
		phi[2]=V2[i][j1][k]  -V2[i1][j][k1];
		phi[3]=V2[i][j1][k1] -V2[i1][j][k];
		dv21=(phi[0]+phi[1]+phi[2]+phi[3])/dl;
		dv22=(phi[0]+phi[1]-phi[2]-phi[3])/dl;
		dv23=(phi[0]-phi[1]+phi[2]-phi[3])/dl;


		phi[0]=V3[i][j][k]   -V3[i1][j1][k1];
		phi[1]=V3[i][j][k1]  -V3[i1][j1][k];
		phi[2]=V3[i][j1][k]  -V3[i1][j][k1];
		phi[3]=V3[i][j1][k1] -V3[i1][j][k];
		dv31=(phi[0]+phi[1]+phi[2]+phi[3])/dl;
		dv32=(phi[0]+phi[1]-phi[2]-phi[3])/dl;
		dv33=(phi[0]-phi[1]+phi[2]-phi[3])/dl;

		for(I=0;I<6;I++){
			ds[I]=
			 Cij[I][0]*dv11+Cij[I][1]*dv22+Cij[I][2]*dv33
			+Cij[I][3]*(dv32+dv23)
			+Cij[I][4]*(dv13+dv31)
			+Cij[I][5]*(dv21+dv12);
		}
		S1[i][j][k]+=ds[0];
		S2[i][j][k]+=ds[1];
		S3[i][j][k]+=ds[2];
		S4[i][j][k]+=ds[3];
		S5[i][j][k]+=ds[4];
		S6[i][j][k]+=ds[5];
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

	double dl=4.*dh/dt*rho;

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
		ds11=(phi[0]+phi[1]+phi[2]+phi[3])/dl;

		phi[0]=S2[i1][j1][k1]-S2[i][j][k];
		phi[1]=S2[i1][j1][k] -S2[i][j][k1];
		phi[2]=S2[i1][j][k1] -S2[i][j1][k];
		phi[3]=S2[i1][j][k]  -S2[i][j1][k1];
		ds22=(phi[0]+phi[1]-phi[2]-phi[3])/dl;

		phi[0]=S3[i1][j1][k1]-S3[i][j][k];
		phi[1]=S3[i1][j1][k] -S3[i][j][k1];
		phi[2]=S3[i1][j][k1] -S3[i][j1][k];
		phi[3]=S3[i1][j][k]  -S3[i][j1][k1];
		ds33=(phi[0]-phi[1]+phi[2]-phi[3])/dl;
	//	if(abs(ds33)>0.0) printf("ds33=%lf\n",ds33);
	//	if(abs(S3[i1][j1][k1])>0.0) printf("S3=%lf,phi=%lf %lf %lf %lf dl=%lf\n",S3[i1][j1][k1],phi[0],phi[1],phi[2],phi[3],dl);

		phi[0]=S4[i1][j1][k1]-S4[i][j][k];
		phi[1]=S4[i1][j1][k] -S4[i][j][k1];
		phi[2]=S4[i1][j][k1] -S4[i][j1][k];
		phi[3]=S4[i1][j][k]  -S4[i][j1][k1];
		ds42=(phi[0]+phi[1]-phi[2]-phi[3])/dl;
		ds43=(phi[0]-phi[1]+phi[2]-phi[3])/dl;

		phi[0]=S5[i1][j1][k1]-S5[i][j][k];
		phi[1]=S5[i1][j1][k] -S5[i][j][k1];
		phi[2]=S5[i1][j][k1] -S5[i][j1][k];
		phi[3]=S5[i1][j][k]  -S5[i][j1][k1];
		ds51=(phi[0]+phi[1]+phi[2]+phi[3])/dl;
		ds53=(phi[0]-phi[1]+phi[2]-phi[3])/dl;

		phi[0]=S6[i1][j1][k1]-S6[i][j][k];
		phi[1]=S6[i1][j1][k] -S6[i][j][k1];
		phi[2]=S6[i1][j][k1] -S6[i][j1][k];
		phi[3]=S6[i1][j][k]  -S6[i][j1][k1];
		ds61=(phi[0]+phi[1]+phi[2]+phi[3])/dl;
		ds62=(phi[0]+phi[1]-phi[2]-phi[3])/dl;

		V1[i][j][k]+=(ds11+ds62+ds53);
		V2[i][j][k]+=(ds61+ds22+ds43);
		V3[i][j][k]+=(ds51+ds42+ds33);
	}
	}
	}
};

void DOMAIN::set_size(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];
	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf %lf\n",Xa,Xa+1,Xa+2);
	fscanf(fp,"%lf %lf %lf\n",Xb,Xb+1,Xb+2);
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&n_lmb);
	fclose(fp);
	printf("Xb=%lf %lf %lf\n",Xb[0],Xb[1],Xb[2]);
	printf("nlmb=%d\n",n_lmb);
	fclose(fp);
};
void DOMAIN::set_cij(char *fname){
	src.load_prms(fname);
	cij.fload(fname);
	cij.print_cij();
	rho=cij.rho;
	printf("rho=%lf\n",rho);
};
void DOMAIN::set_src(char *fname){
	src.load_prms(fname);

	int ndim=3;
	int n_T=1.5*ceil( cij.cL/cij.cT*n_lmb*sqrt(ndim));
	if( n_T > src.n_T) src.n_T=n_T;

	double CFL;

	dh=src.T0*cij.cT/n_lmb;
	dt=src.T0/src.n_T;
	CFL=cij.cL*dt*sqrt(ndim)/dh;
	printf("CFL=%lf\n",CFL);
	printf("dt=%lf, dh=%lf\n",dt,dh);
	printf("n_T=%d, n_lmb=%d\n",src.n_T,n_lmb);

	src.set_wvfm();

};


void DOMAIN::setup(
){
	int i;
	for(i=0;i<3;i++){
		Wd[i]=Xb[i]-Xa[i];
		Ndiv[i]=ceil(Wd[i]/dh);
		Wd[i]=dh*Ndiv[i];
		Xb[i]=Xa[i]+dh*Ndiv[i];
	};

	printf("Xa=%lf %lf %lf\n",Xa[0],Xa[1],Xa[2]);
	printf("Xb=%lf %lf %lf\n",Xb[0],Xb[1],Xb[2]);
	printf("Ndiv=%d %d %d\n",Ndiv[0],Ndiv[1],Ndiv[2]);

	rho=cij.rho;
	printf("rho=%lf\n",rho);

	fld.mem_alloc(Ndiv);
	fld.dh=dh;
	fld.rho=rho;
	fld.Cij=cij.cij;

	dt=src.wv.dt;
	Nt=src.wv.Nt;
	fld.dt=dt;
	printf("CFL=%lf\n",Courant(dt,cij.cL,dh));


	int j,k,nsg=0,Ng[3];
	int indx1[3],indx2[3];
	cod2indx(src.xs1,indx1,0);
	cod2indx(src.xs2,indx2,0);
	Ng[0]=Ndiv[0]+1;
	Ng[1]=Ndiv[1]+1;
	Ng[2]=Ndiv[2]+1;
	for(i=0;i<3;i++){
		if(indx1[i]<0) indx1[i]=0;
		if(indx2[i]>=Ng[i]-1) indx1[i]=Ng[i]-1;
	};

	for(i=0;i<3;i++){
		src.indx1[i]=indx1[i];
		src.indx2[i]=indx2[i];
	};
	for(i=indx1[0];i<indx2[0]+1;i++){
	for(j=indx1[1];j<indx2[1]+1;j++){
	for(k=indx1[2];k<indx2[2]+1;k++){
		nsg++;
	}
	}
	}
	printf("indx1=%d %d %d\n",indx1[0],indx1[1],indx1[2]);
	printf("indx2=%d %d %d\n",indx2[0],indx2[1],indx2[2]);
	printf("nsg=%d\n",nsg);
		
};
void DOMAIN::cod2indx(
	double *xcod,	// coordinate (x,y)
	int *indx,	// index (i,j)
	int type	// 0:stress, 1: velocity
){
	if(type==0){
		indx[0]=floor((xcod[0]-Xa[0])/dh);
		indx[1]=floor((xcod[1]-Xa[1])/dh);
		indx[2]=floor((xcod[2]-Xa[2])/dh);
	}else if(type==1){
		indx[0]=floor((xcod[0]-Xa[0])/dh+0.5);
		indx[1]=floor((xcod[1]-Xa[1])/dh+0.5);
		indx[2]=floor((xcod[2]-Xa[2])/dh+0.5);
	}else{
		puts("Invalid grid type specified in cod2indx");
		exit(-1);
	}
};
double Courant(double dt, double vel, double ds){
	return(vel*dt/ds*sqrt(3.0));
};
void DOMAIN::apply_source(int it){
	int i,j,k;
	int indx1[3],indx2[3];
	for(i=0;i<3;i++){
		indx1[i]=src.indx1[i];
		indx2[i]=src.indx2[i];
	};
	double amp=src.wv.amp[it];
	for(i=indx1[0]; i<=indx2[0]; i++){
	for(j=indx1[1]; j<=indx2[1]; j++){
	for(k=indx1[2]; k<=indx2[2]; k++){
		fld.S3[i][j][k]=amp;
	}
	}
	}
};

void DOMAIN::write_xslice(int it,double xout){
	static int nout=0;
	char fname[128];
	sprintf(fname,"v%dx.out",nout);
	FILE *fp=fopen(fname,"w");

	double v1,v2,v3;
	int i,j,k;
	i=int((xout-Xa[0])/dh);

	while(i<0) i+=Ndiv[0];
	while(i>=Ndiv[0]) i-=Ndiv[0];
	xout=(i+0.5)*dh+Xa[0];

	fprintf(fp,"# time\n");
	fprintf(fp,"%lf\n",it*dt);
	fprintf(fp,"# Xa[0:3]\n");
	fprintf(fp,"%lf, %lf, %lf\n",xout,Xa[1],Xa[2]);
	fprintf(fp,"# Xb[0:3]\n");
	fprintf(fp,"%lf, %lf, %lf\n",xout,Xb[1],Xb[2]);
	fprintf(fp,"# Ndiv[0:3]\n");
	fprintf(fp,"%d, %d, %d\n",1,Ndiv[1],Ndiv[2]);
	fprintf(fp,"# v1, v2, v3 \n");

	for(j=0;j<Ndiv[1];j++){
	for(k=0;k<Ndiv[2];k++){
		v1=fld.V1[i][j][k];
		v2=fld.V2[i][j][k];
		v3=fld.V3[i][j][k];
		fprintf(fp,"%lf, %lf, %lf\n",v1,v2,v3);
	}
	}
	fclose(fp);
	nout++;
};
void DOMAIN::write_zslice(int it,double zout){
	static int nout=0;
	char fname[128];
	sprintf(fname,"v%dz.out",nout);
	FILE *fp=fopen(fname,"w");

	double v1,v2,v3;
	int i,j,k;
	k=int((zout-Xa[2])/dh);

	while(k<0) k+=Ndiv[2];
	while(k>=Ndiv[2]) k-=Ndiv[2];
	zout=(k+0.5)*dh+Xa[2];

	fprintf(fp,"# time\n");
	fprintf(fp,"%lf\n",it*dt);
	fprintf(fp,"# Xa[0:3]\n");
	fprintf(fp,"%lf, %lf, %lf\n",Xa[0],Xa[1],zout);
	fprintf(fp,"# Xb[0:3]\n");
	fprintf(fp,"%lf, %lf, %lf\n",Xb[0],Xb[1],zout);
	fprintf(fp,"# Ndiv[0:3]\n");
	fprintf(fp,"%d, %d, %d\n",Ndiv[0],Ndiv[1],1);
	fprintf(fp,"# v1, v2, v3 \n");

	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		v1=fld.V1[i][j][k];
		v2=fld.V2[i][j][k];
		v3=fld.V3[i][j][k];
		fprintf(fp,"%lf, %lf, %lf\n",v1,v2,v3);
	}
	}
	fclose(fp);
	nout++;
};
void DOMAIN::write_yslice(int it,double yout){
	static int nout=0;
	char fname[128];
	sprintf(fname,"v%dy.out",nout);
	FILE *fp=fopen(fname,"w");

	double v1,v2,v3;
	int i,j,k;
	j=int((yout-Xa[1])/dh);

	while(j<0) j+=Ndiv[1];
	while(j>=Ndiv[1]) j-=Ndiv[1];
	yout=(j+0.5)*dh+Xa[1];

	fprintf(fp,"# time\n");
	fprintf(fp,"%lf\n",it*dt);
	fprintf(fp,"# Xa[0:3]\n");
	fprintf(fp,"%lf, %lf, %lf\n",Xa[0],yout,Xa[2]);
	fprintf(fp,"# Xb[0:3]\n");
	fprintf(fp,"%lf, %lf, %lf\n",Xb[0],yout,Xb[2]);
	fprintf(fp,"# Ndiv[0:3]\n");
	fprintf(fp,"%d, %d, %d\n",Ndiv[0],1,Ndiv[2]);
	fprintf(fp,"# v1, v2, v3 \n");

	for(i=0;i<Ndiv[0];i++){
	for(k=0;k<Ndiv[2];k++){
		v1=fld.V1[i][j][k];
		v2=fld.V2[i][j][k];
		v3=fld.V3[i][j][k];
		fprintf(fp,"%lf, %lf, %lf\n",v1,v2,v3);
	}
	}
	fclose(fp);
	nout++;
};

void DOMAIN::write_v(char *fname, int it){
	FILE *fp=fopen(fname,"w");
	int i,j,k;
	k=Ndiv[2]-1;
	double v1,v2,v3;
	fprintf(fp,"# time\n");
	fprintf(fp,"%lf\n",it*dt);
	fprintf(fp,"# Xa[0:3]\n");
	fprintf(fp,"%lf, %lf, %lf\n",Xa[0],Xa[1],Xa[2]);
	fprintf(fp,"# Xb[0:3]\n");
	fprintf(fp,"%lf, %lf, %lf\n",Xb[0],Xb[1],Xb[2]);
	fprintf(fp,"# Ndiv[0:3]\n");
	fprintf(fp,"%d, %d, %d\n",Ndiv[0],Ndiv[1],Ndiv[2]);
	fprintf(fp,"# v1, v2, v3 \n");
	/*
	for(i=0;i<Ndiv[0];i++){
	for(j=0;j<Ndiv[1];j++){
		v1=fld.V1[i][j][k];
		v2=fld.V2[i][j][k];
		v3=fld.V3[i][j][k];
		fprintf(fp,"%lf, %lf, %lf\n",v1,v2,v3);
	}
	}
	*/
	j=Ndiv[1]/2;
	for(i=0;i<Ndiv[0];i++){
	for(k=0;k<Ndiv[2];k++){
		v1=fld.V1[i][j][k];
		v2=fld.V2[i][j][k];
		v3=fld.V3[i][j][k];
		fprintf(fp,"%lf, %lf, %lf\n",v1,v2,v3);
	}
	}
	fclose(fp);
};
