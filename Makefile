COfit.o: COfit.cpp COfit.h lib/idlwhere.h
	gcc -fpermissive -std=c++11 -o COfit  COfit.cpp -Ofast  -larmadillo -lstdc++ -lm -llapack -lblas

