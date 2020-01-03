class STIFF{
	public:
		int type;
//		double cij[6][6];
		double **cij;
		double lmb,mu;
		double cL,cT,rho;
		void load(int i);
		void fload(char *fname);
		STIFF();
		void print_cij();
	private:
};

