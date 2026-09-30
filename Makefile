CC ?= gcc
CFLAGS ?= -std=c11 -Wall -Wextra -O2
SRCS = field.c link.c source.c relay.c dest.c sim.c main.c

.PHONY: all run clean

all: bats_sim

bats_sim: $(SRCS) bats.h
	$(CC) $(CFLAGS) -o bats_sim $(SRCS)

run: bats_sim
	./bats_sim input.txt recovered.txt

clean:
	rm -f bats_sim bats_sim.exe recovered.txt
