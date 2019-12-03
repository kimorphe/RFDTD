class Cmplx{	// Complex Number Class
	public:
		double re,im;
		double Abs();
		Cmplx(double a, double b);
		Cmplx();
		Cmplx Cnjg();
		void disp();
		void set(double a, double b);
		Cmplx operator+(Cmplx z);
		Cmplx operator-(Cmplx z);
		Cmplx operator*(Cmplx z);
		Cmplx operator/(Cmplx z);
		Cmplx operator*(double z);
		void operator=(double z);
	private:
};
Cmplx exp(Cmplx z);
class InWv{
	public :
		double t1,t2,dt,T0;
		double *amp;
		int Nt,wvtyp,nbrst,nwv;
		int nsig; //narrowness factor used in Gaussian amplitude modulation
		InWv(char *);	// constructor (data read from file char*) 
		InWv(int Nt);	// 
		void disp();
		void gen_sin();		
		void gen_cos();		
		void gen_wvfm();
		void set_taxis(double t1, double t2);
		void out(char *);
		double df,fmax;
		Cmplx *Amp;
		void DFT();
		void DFTout(char *);
		void Amod(double tb, int nsig);
	private:
		void mem_alloc();
};		
