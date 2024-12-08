cc := /user/bin/gcc
RM := /usr/bin/rm

CFLAGS := -c -std=c99 -Wall -Wextra -ggdb
LFLAGS:= -std=c99 -O3 -ggdb

SRC := $(wildcard *.c)
OBJS:= $(SRC:.c=.o)

.PHONY : all clean
all: MD

MD: $(OBJS)
	gcc -fopenmp $(LFLAGS) $(OBJS) -o $@ -lm
%.o: %.c
	gcc -fopenmp $(CFLAGS) -c $< -o $@ -I.

clean:
	$(RM) -f *.o MD