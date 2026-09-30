CC ?= gcc
CFLAGS ?= -std=c11 -Wall -Wextra -O2
SRCS = field.c link.c source.c relay.c dest.c codec.c sim.c main.c
LIB_SRCS = field.c source.c dest.c codec.c

.PHONY: all lib run clean

all: bats_sim

bats_sim: $(SRCS) bats.h bats_codec.h
	$(CC) $(CFLAGS) -o bats_sim $(SRCS)

lib: libbats.a

libbats.a: $(LIB_SRCS) bats.h bats_codec.h
	$(CC) $(CFLAGS) -c $(LIB_SRCS)
	ar rcs libbats.a field.o source.o dest.o codec.o

run: bats_sim
	./bats_sim input.txt recovered.txt

clean:
	rm -f bats_sim bats_sim.exe recovered.txt libbats.a *.o
