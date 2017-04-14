filter.svg: main Makefile lpf-140MHz.txt
	gnuplot filter.gnuplot
lpf-140MHz.txt: main
	./main >lpf-140MHz.txt
main: main.cpp Makefile
	g++ -o main main.cpp --std=c++11
