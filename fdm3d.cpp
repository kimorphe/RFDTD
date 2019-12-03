#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cij.h"
#include "field.h"

#define __field 1

int main(int argv, char *argc[]){
	STIFF cij;

	cij.load(0);
	cij.print_cij();

	Field fld;
	int Ndiv[3]={2,3,4};
	fld.mem_alloc(Ndiv);
	fld.print_Vt();

	DOMAIN dom;
	double Xa[3]={0.0,0.0,0.0};
	double Xb[3]={20.0,30.0,10.0};


	dom.setup(Xa,Xb,0.1);

	return(0);
};
