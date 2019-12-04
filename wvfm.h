#include <complex>
using namespace std;
class InWv{
	public :
		double t1,t2,dt,T0;
		double *amp;
		int Nt,wvtyp,nbrst,nwv;
		int nsig; //narrowness factor used in Gaussian amplitude modulation
		void read_prms(char *);	// constructor (data read from file char*) 
		void set_Nt(int Nt);	// 
		InWv(int Nt);	// 
		InWv();	// 
		void disp();
		void gen_sin();		
		void gen_cos();		
		void gen_wvfm();
		void set_taxis(double t1, double t2);
		void out(char *);
		double df,fmax;
		complex<double> *Amp;
		void DFT();
		void DFTout(char *);
		void Amod(double tb, int nsig);
	private:
		void mem_alloc();
};		
