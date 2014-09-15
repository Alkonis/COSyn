COfit.o: COfit.cpp COfit.h lib/idlarma.h include/armadillo
	mpicc -fpermissive -std=c++11 -o COfit  COfit.cpp -Ofast -lstdc++ -L/usr/lib/gcc/x86_64-linux-gnu/4.7 -L./include -lopenblas  -lgfortran -lm -pthread 

