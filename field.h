//#include"wvfm.h"
#include"source.h"
double Courant(double dt, double vel, double ds);
class Field{
	public:
		double ***V1,***V2,***V3;
		double ***S1,***S2,***S3;
		double ***S4,***S5,***S6;
		void init();
		double ***mem_alloc3d(int Ndiv[3]);
		void mem_alloc(int Ndiv[3]);
		void print_Vt();
		void s2v();
		void v2s();
		int Ndiv[3];
		int ngs[3],ngv[3];
		double dh,dt,rho;
		//double Cij[6][6];
		double **Cij;
	private:
};
class DOMAIN{
	public:
		double Xa[3],Xb[3];
		int Ndiv[3];
		int Nt;
		double dh,dt,dx[3];
		double rho;
		STIFF cij;
		Field fld;
		void setup(double Xa[3],double Xb[3],double dh);
//		InWv wv;
		SRC src;
		void cod2indx(double *xcod, int *indx, int type);
		void apply_source(int it);
		void write_v(char *fname);
	private:
};

