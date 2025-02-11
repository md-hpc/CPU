
gxx = g++ -lm common.cpp
common = common.cpp common.h
bin = reference list cells hemisphere hemilist

all: stamp-all build

debug: gxx += -g -DDEBUG
debug: stamp-debug build

stamp-debug:
	-rm stamp-all $(bin)
	touch stamp-debug
stamp-all:
	-rm stamp-debug $(bin)
	touch stamp-all

clean:
	-rm $(bin) test

build: $(bin)

%: %.cpp $(common)
	$(gxx) $@.cpp -o $@
