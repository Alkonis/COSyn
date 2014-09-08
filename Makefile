COfit.o: COfit.cpp COfit.h lib/idlwhere.h include/armadillo
	gcc -static -fpermissive -std=c++11 -o COfit  COfit.cpp -Ofast -lstdc++ -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.7 -llapack -lgfortran -lopenblas -lm -larpack -pthread 


