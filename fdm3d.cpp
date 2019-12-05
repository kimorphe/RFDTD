#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cij.h"
#include "field.h"

#define __field 1

int main(int argv, char *argc[]){

	DOMAIN dom;
	double Xa[3]={0.0,0.0,0.0};
	double Xb[3]={20.0,30.0,10.0};
	double dx=0.1;
	dom.setup(Xa,Xb,dx);

	int i,num=1;
	char fname[128];
	FILE *fp;
	for(i=0;i<dom.Nt;i++){
		printf("i=%d\n",i);
		dom.apply_source(i);
		dom.fld.s2v();
		dom.fld.v2s();

		sprintf(fname,"v%d.out",num);
		puts(fname);
		dom.write_v(fname);
		num++;
	};

	return(0);
};

