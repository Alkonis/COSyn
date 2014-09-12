COfit.o: COfit.cpp COfit.h lib/idlarma.h include/armadillo
	gcc -fpermissive -std=c++11 -o COfit  COfit.cpp -Ofast -lstdc++ -Bstatic -llapack -L/usr/lib/gcc/x86_64-linux-gnu/4.7 -lgfortran -lopenblas -larpack -lm -pthread -Bdynamic

