fdm3d:	fdm3d.o cij.o field.o wvfm.o source.o recs.o
	g++ -o fdm3d fdm3d.o cij.o field.o wvfm.o source.o recs.o

fdm3d.o: fdm3d.cpp 
	g++ -c fdm3d.cpp
field.o: field.cpp field.h
	g++ -c field.cpp
cij.o:	cij.cpp cij.h
	g++ -c cij.cpp

wvfm.o: wvfm.cpp wvfm.h
	g++ -c wvfm.cpp

source.o: source.cpp source.h
	g++ -c source.cpp
recs.o: recs.cpp recs.h
	g++ -c resc.cpp
