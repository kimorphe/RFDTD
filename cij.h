class STIFF{
	public:
		int type;
		double cij[6][6];
		double lmb,mu;
		void load(int i);
		STIFF();
		void print_cij();
	private:
};

