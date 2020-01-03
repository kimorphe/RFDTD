#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cij.h"

#define DB 0

STIFF::STIFF(){
	int i,j;
	double *pt=(double *)malloc(sizeof(double)*6*6);
	cij=(double **)malloc(sizeof(double *)*6);

	for(i=0;i<6;i++){
		cij[i]=pt+i*6;
	};
	

	for(i=0;i<6;i++){
	for(j=0;j<6;j++){
		cij[i][j]=0.0;
	}
	}
};
void STIFF::fload(char *fname){
	char cbff[128];
	FILE *fp=fopen(fname,"r");
	if(fp==NULL){
		printf("Cannot find file %s\n",fname);	
		printf(" --> abort process.\n");
		exit(-1);
	}
	
	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&type);

	if(type==0){
		fgets(cbff,128,fp);
		fscanf(fp,"%lf\n",&rho);
		fgets(cbff,128,fp);
		fscanf(fp,"%lf %lf\n",&cL,&cT);

		mu=rho*cT*cT;	// [GPa]
		lmb=rho*(cL*cL-2.*cT*cT); //[GPa]
		double K=lmb+2.*mu;
		cij[0][0]=K;
		cij[1][1]=K;
		cij[2][2]=K;
		cij[3][3]=mu;
		cij[4][4]=mu;
		cij[5][5]=mu;

		cij[0][1]=lmb; cij[1][0]=lmb;
		cij[0][2]=lmb; cij[2][0]=lmb;
		
	}else{
		printf("Thie version cannot read anisotropic stiffness tensor components\n");
		exit(-1);
	}
	
	fclose(fp);
};
void STIFF::load(int type){
	if(type==0){	// isotropic material
		cL=6.0;	//[km/s]=[micro sec/mm]
		cT=3.0;	//[km/s]=[micro sec/mm]
		rho=7.9; //[g/cm^3]=[1,000 x kg/m^3]
		mu=rho*cT*cT;	// [GPa]
		lmb=rho*(cL*cL-2.*cT*cT); //[GPa]
		double K=lmb+2.*mu;
		cij[0][0]=K;
		cij[1][1]=K;
		cij[2][2]=K;
		cij[3][3]=mu;
		cij[4][4]=mu;
		cij[5][5]=mu;

		cij[0][1]=lmb; cij[1][0]=lmb;
		cij[0][2]=lmb; cij[2][0]=lmb;
		cij[1][2]=lmb; cij[2][1]=lmb;
	};


};

void STIFF::print_cij(){
	int i,j;
	printf("cij: stiffness constant in GPa\n[\n");
	for(i=0;i<6;i++){
	for(j=0;j<6;j++){
		printf("%lf ",cij[i][j]);
	}
	printf("\n");
	}
	printf("]\n");
};

#if DB == 1
int main(int argv, char *argc[]){
	STIFF cij;

	cij.load(0);
	cij.print_cij();

	return(0);
};
#endif
