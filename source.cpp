#include "source.h"

void SRC::setup(){

	int Nt=300;
	wv.set_Nt(Nt);
	wv.set_taxis(0.0,1.0);	// [t1,t2]
	wv.T0=0.2;		// T0 (period)
	wv.gen_wvfm();		// generate waveform
	wv.Amod(0.1,4);		// amlitude modlulation

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
