
gxx = g++ -lm common.cpp
common = common.cpp common.h
bin = reference list cells hemisphere

all: build

debug_cxxflags: gxx += -g -DDEBUG

debug: debug_cxxflags build

clean:
	-rm $(bin)

build: $(bin)


%: %.cpp $(common)
	$(gxx) $@.cpp -o $@
