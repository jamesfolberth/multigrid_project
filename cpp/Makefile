# GNU Makefile

# use normal g++ if you don't have GCC version >= 4.9
#CPP=g++ 
CPP=g++ -fdiagnostics-color=auto
SHELL=/usr/bin/bash

#CPPFLAGS=-Wall -Wpedantic -std=c++11 -g -march=native
CPPFLAGS=-Wall -Wpedantic -std=c++11 -pg -fprofile-arcs -ftest-coverage -Og -march=native
#CPPFLAGS=-Wall -Wpedantic -std=c++11 -O2 -march=native

LDFLAGS=
LDLIBS=-lc -lm

OPENCVLIBS=`pkg-config --cflags --libs opencv`

COMMON_CPP= matrix_coo.cpp\
				matrix_crs.cpp\
				utils.cpp\
				model_problems.o\
				classical_solvers.cpp\
				multigrid.cpp

COMMON_HPP= $(COMMON_CPP:%.cpp=%.hpp)
COMMON_O  = $(COMMON_CPP:%.cpp=%.o)

.PHONY: all clean run

all: test run

run: test
	./test

# Test
test: $(COMMON_O) test.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c test.cpp
	$(CPP) $(CPPFLAGS) -flto $(LDLIBS) $(LDFLAGS) $(COMMON_O) test.o -o test


# ./testing
display_image: testing/display_image.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(OPENCVLIBS) $(LDFLAGS) testing/display_image.cpp -o display_image

# Common files
matrix_coo.o: matrix_coo.cpp matrix_coo.hpp matrix_crs.cpp matrix_crs.hpp utils.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c matrix_coo.cpp

matrix_crs.o: matrix_coo.cpp matrix_coo.hpp matrix_crs.cpp matrix_crs.hpp utils.cpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c matrix_crs.cpp

utils.o: utils.cpp utils.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c utils.cpp

model_problems.o: matrix_coo.cpp matrix_crs.cpp utils.cpp model_problems.cpp model_problems.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c model_problems.cpp

classical_solvers.o: matrix_coo.cpp matrix_crs.cpp utils.cpp classical_solvers.cpp classical_solvers.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c classical_solvers.cpp

multigrid.o: matrix_crs.cpp utils.cpp classical_solvers.cpp multigrid.cpp multigrid.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c multigrid.cpp


# Image Segmentation
image_seg.o: $(COMMON_O) image_seg.cpp image_seg.hpp
	$(CPP) $(CPPFLAGS) $(LDLIBS) $(LDFLAGS) -c image_seg.cpp

image_seg: $(COMMON_O) image_seg.o
	$(CPP) $(CPPFLAGS) -flto $(LDLIBS) $(OPENCVLIBS) $(LDFLAGS) $(COMMON_O) image_seg.o -o image_seg


force:

clean:
	@rm -f *.o
	@rm -f test display_image image_seg
	@rm -f *.out # gprof

cleandata:
	make -C gen_imgs/ cleandata

deepclean: clean cleandata
	#make -C figures/ deepclean

