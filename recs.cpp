#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"recs.h"


/*
class REC{
	public:
		double xcod[3];	// coordinate
		int indx[3];	// grid index
		double **v;	// velocity
		int type;	// 0:stress, 1: velocity
		int Nt;		// time steps
		void set_indx(double Xa[3], double dx);
		void mem_alloc(int Nt);
	private:
};
*/
void REC::set_indx(double Xa[3], double dh){
	if(type==0){	// Stress grid
		indx[0]=floor((xcod[0]-Xa[0])/dh);
		indx[1]=floor((xcod[1]-Xa[1])/dh);
		indx[2]=floor((xcod[2]-Xa[2])/dh);
	}else if(type==1){	// Velocity grid
		indx[0]=floor((xcod[0]-Xa[0])/dh+0.5);
		indx[1]=floor((xcod[1]-Xa[1])/dh+0.5);
		indx[2]=floor((xcod[2]-Xa[2])/dh+0.5);
	}else{
		puts("Invalid grid type specified in cod2indx");
		exit(-1);
	}
};
void REC::mem_alloc(int n){
	Nt=n;
	double *pt=(double *)malloc(sizeof(double)*3*Nt);
	v=(double **)malloc(sizeof(double *)*3);
	for(int i=0;i<3*Nt;i++) pt[i]=0.0;
	for(int i=0;i<3;i++) v[i]=pt+i*Nt;
};
