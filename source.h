#include "wvfm.h"
class SRC{
	public:
		double	xs1[3],xs2[3];	// source aperture (rectangle)
		int indx1[3],indx2[3];
		int *ksrc;	// list of source grid No.s
		int nsg;	// number of source grid
		int stype[6];	// source type (0:off, 1: on)
		InWv wv;	// waveform (source function)
		void setup();
	private:
};
