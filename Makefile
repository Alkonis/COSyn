COfit.o: COfit.cpp COfit.h lib/idlwhere.h
	gcc -fpermissive -std=c++11 -o COfit  COfit.cpp -O2 -larmadillo -lstdc++ -lm -llapack -lblas

