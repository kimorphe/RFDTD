#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "wvfm.h"
#include <complex.h>

#define DB 0

using namespace std;

//-----------------------------------------------------------
InWv::InWv(){
};
InWv::InWv(int N){	// empty constructor
	t1=0.0; t2=1.0;	// time range
	dt=0.0; // sampling interval
	T0=0.1;	// fundamental period
	if(N>1) dt=(t2-t1)/(N-1);
	wvtyp=1; 	// 1:sine, 2: cosine
	nbrst=3;	// burst cycles
	Nt=N;
	mem_alloc();	
	
};
void InWv::set_Nt(int N){	// empty constructor
	t1=0.0; t2=1.0;	// time range
	dt=0.0; // sampling interval
	T0=0.1;	// fundamental period
	if(N>1) dt=(t2-t1)/(N-1);
	wvtyp=1; 	// 1:sine, 2: cosine
	nbrst=3;	// burst cycles
	Nt=N;
	mem_alloc();	
}
	
void InWv::read_prms(char *fname){
	int i;
	FILE *fp;
	char ch[128];
	fp=fopen(fname,"r");
	fgets(ch,128,fp);
	fscanf(fp,"%d\n",&wvtyp);	

	fgets(ch,128,fp);
	fscanf(fp,"%d\n",&Nt);	

	fgets(ch,128,fp);
	fscanf(fp,"%lf %d\n",&T0,&nbrst);	

	fgets(ch,128,fp);
	fscanf(fp,"%d\n",&nsig); // narrowness factor used in Gaussian amp. modulation	

	fgets(ch,128,fp);
	fscanf(fp,"%lf %lf\n",&t1,&t2);	

	dt=(t2-t1)/(Nt-1);

	mem_alloc();
	gen_wvfm();

	double tb=nbrst*T0*0.5;
	if(nsig != 0) Amod(tb,nsig);
	if(wvtyp <0){	// read numerical waveform if wvtyp is negative
		fgets(ch,128,fp);
		for(i=0;i<Nt;i++){
			 fscanf(fp,"%lf\n",amp+i);	
		}
	}
	fclose(fp);
};

void InWv :: disp(){
	cout<<"**** waveform parameters ****\n";
	cout<<"wtyp="<<wvtyp<<'\n';
	cout<<"t1="<<t1<<" "<<"t2="<<t2<<'\n';
	cout<<"dt="<<dt<<'\n';
	cout<<"Nt="<<Nt<<'\n';
	cout<<"nbrst="<<nbrst<<'\n';
	cout<<"T0="<<T0<<'\n';
}

void InWv::mem_alloc(){
	amp=(double *)malloc(sizeof(double)*Nt);
}

void InWv::set_taxis(double s1, double s2){
	t1=s1;
	t2=s2;
	dt=(t2-t1)/int(Nt-1);
}

void InWv:: gen_wvfm(){
	switch(wvtyp){
	case 1: // sine function
		gen_sin();
		break;
	case 2:	// cosine function
		gen_cos();
		break;
	};
}

void InWv:: gen_sin(){
	double PI=4.0*atan(1.0);
	double omg=2.*PI/T0,tt;

	for(int i=0;i<Nt;i++){
		tt=t1+dt*i;
		amp[i]=0.0;
		if(tt-t1<=nbrst*T0) amp[i]=sin(omg*tt);
	}
}
void InWv:: gen_cos(){
	double PI=4.0*atan(1.0);
	double omg=2.*PI/T0,tt;

	for(int i=0;i<Nt;i++){
		tt=t1+dt*i;
		amp[i]=0.0;
//		if(tt-t1<=nbrst*T0) amp[i]=1-cos(omg*tt);
		if(tt-t1<=nbrst*T0) amp[i]=cos(omg*tt);
	}
}
void InWv::out(char *fout){
	FILE *fp=fopen(fout,"w");
	int i;
	for(i=0;i<Nt;i++) fprintf(fp,"%lf %lf\n",t1+dt*i, amp[i]);
	fclose(fp);
}


void InWv::Amod(
	double tb,	// mean 
	int nsig)	// narrowness 
{
	int i;
	double arg,sig=tb/nsig;
	for(i=0;i<Nt;i++){
		arg=(i*dt+t1-tb)/sig;
		arg*=arg;
		amp[i]*=exp(-arg*0.5);
	
	}
};
//---------------------------------------------------

void InWv::DFT(){
	int i,j;
	double omg,pi=4.0*atan(1.0),tt;
	fmax=1./dt;
	df=fmax/Nt;	
	complex <double>zi(0.0,1.0);
	double dw=2.0*pi*df;

	Amp=(complex<double> *)malloc(sizeof(complex<double>)*Nt);
	for(i=0;i<Nt;i++){
		Amp[i]=complex<double>(0.0,0.0);
		tt=t1+dt*i;
	for(j=0;j<Nt;j++){
		omg=dw*j;
		Amp[i]=Amp[i]+exp(zi*omg*tt)*amp[j];
	}	
	}

}
void InWv::DFTout(char *fout){
	FILE *fp=fopen(fout,"w");
	int i;
	for(i=0;i<Nt;i++) fprintf(fp,"%lf %lf %lf %lf \n",df*i, Amp[i].real(),Amp[i].imag(),abs(Amp[i]));
	fclose(fp);
}

#if DB == 1
int main(){
	char fname[]="wvfm.out";
	char Fname[]="wvfm.dft";
	InWv wv(200);
	wv.set_taxis(0.0,1.0);
	wv.T0=0.2;
	wv.gen_wvfm();
	wv.Amod(0.1,3);
	wv.out(fname);

	wv.DFT();
	wv.DFTout(Fname);

	complex <double> z,zi,w;
	zi=complex<double>(0.0,1.0);
	z=complex<double>(2.0,3.0);
	w=z*zi;
	printf("%lf %lf\n",w.real(),w.imag());


	return(0);
};
#endif
