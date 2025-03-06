
src = common.cpp particle.cpp particle8.cpp workers.cpp simulation.cpp simops.cpp avx.cpp
gxx = g++ -march=native -lm -g

all:
	$(gxx) -o sim -O3 md.cpp $(src)
debug: 
	$(gxx) -o sim -O0 md.cpp $(src)
test:
	$(gxx) -o tst test.cpp $(src)

