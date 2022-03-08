FLAGS := -O3

main: setup cmf
	mpicxx -g --std=c++11 -I. -I${CMF}/include -I${PTL}/include ${FLAGS} main.cc -o pbc -L${CMF}/lib -lcmf -L${PTL}/lib -lPTL

cmf:
	make -C ${CMF} -f makefile

setup:
	mkdir -p output
	mkdir -p series
	mkdir -p checkpoint

clean:
	rm -f pbc
	rm -r output
	rm -r series
