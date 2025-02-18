include ../Makefile

src = $(wildcard *.cpp)
bin = $(src:.cpp=)

all: $(bin)

%: %.cpp
	$(gxx) $< -o $@
