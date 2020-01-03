#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cij.h"
#include "field.h"

#define __field 1

int main(int argv, char *argc[]){

	char finp[128]="fdm3d.inp";
	char fcij[128]="cij.inp";
	char fsrc[128]="source.inp";

	DOMAIN dom;
	dom.set_size(finp);
	dom.set_cij(fcij);
	dom.set_src(fsrc);
	dom.setup();

	exit(-1);


	int i,num=1;
	char fname[128];
	for(i=0;i<dom.Nt;i++){
		printf("i=%d\n",i);
		dom.apply_source(i);
		dom.fld.s2v();
		dom.fld.v2s();
		if(i%5==0){
		sprintf(fname,"v%d.out",num);
		puts(fname);
		dom.write_v(fname);
		num++;
		}
	};

	return(0);
};

