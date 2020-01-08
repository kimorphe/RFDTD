class REC{
	public:
		double xcod[3];	// coordinate
		int indx[3];	// grid index
		double **v;	// velocity
		int type;	// 0:stress, 1: velocity
		int Nt;		// time steps
		void set_indx(double Xa[3], double dx);
		void mem_alloc(int Nt);
	private:
};
