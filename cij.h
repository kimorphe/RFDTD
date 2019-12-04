class STIFF{
	public:
		int type;
//		double cij[6][6];
		double **cij;
		double lmb,mu;
		double cL,cT,rho;
		void load(int i);
		STIFF();
		void print_cij();
	private:
};

