#include "source.h"

void SRC::load_prms(char *fname){

	char cbff[128];
	FILE *fp=fopen(fname,"r");
	if(fp==NULL){
		printf("File %s not found !\n",fname);
		printf("--> abort process..\n");
		exit(-1);
	};

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&T0);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&Td);

	fgets(cbff,128,fp);
	fscanf(fp,"%d\n",&n_T);

	fgets(cbff,128,fp);
	fscanf(fp,"%lf %lf %lf\n",xs1,xs1+1,xs1+2);
	fscanf(fp,"%lf %lf %lf\n",xs2,xs2+1,xs2+2);

	fclose(fp);
};
void SRC::set_wvfm(){

	double dt=T0/n_T;
	int Nt=ceil(Td/dt)+1;
	wv.set_Nt(Nt);
	wv.set_taxis(0.0,Td);	// [t1,t2]
	wv.T0=T0;		// T0 (period)
	wv.gen_wvfm();		// generate waveform
	wv.Amod(2.*T0,4);		// amlitude modlulation
	char tmp[]="wvfm0.out";
	wv.out(tmp);	

	wv.DFT();
	sprintf(tmp,"wvfft.out");
	wv.DFTout(tmp);
	for(int i=0;i<6;i++) stype[i]=0;
	stype[2]=1;
};
void SRC::setup(){

	double Ndt=35;
	double T0=1.0;

	double dtau=T0/Ndt;
	double t_end=10*T0;

	int Nt=ceil(t_end/dtau);
	//int Nt=300;
	wv.set_taxis(0.0,T0*10);	// [t1,t2]
	wv.set_Nt(Nt);
	//wv.T0=1.0;		// T0 (period)
	wv.T0=T0;		// T0 (period)
	wv.gen_wvfm();		// generate waveform
	//wv.Amod(0.1,4);		// amlitude modlulation
	wv.Amod(T0*0.5,4);		// amlitude modlulation

	char fname[]="wvfm.out";
	wv.out(fname);	

	int i;
	for(i=0;i<6;i++) stype[i]=0;
	stype[2]=1;

	xs1[0]=10.0-0.05;
	xs1[1]=15.0-0.05;
	xs1[2]=10.0;

	xs2[0]=10.0+0.05;
	xs2[1]=15.0+0.05;
	xs2[2]=10.0;


};
