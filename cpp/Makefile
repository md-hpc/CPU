
gxx = g++ -lm common.cpp
common = common.cpp common.h
bin = reference list cells hemisphere

all: build

debug: gxx += -g -DDEBUG
debug: build

clean:
	-rm $(bin) test

build: $(bin)


%: %.cpp $(common)
	$(gxx) $@.cpp -o $@
