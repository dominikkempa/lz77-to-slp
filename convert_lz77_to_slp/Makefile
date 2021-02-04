SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -O3 -DNDEBUG -march=native -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2

all: convert_lz77_to_slp

convert_lz77_to_slp:
	$(CC) $(CFLAGS) -o convert_lz77_to_slp main.cpp src/utils.cpp src/karp_rabin_hashing.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f convert_lz77_to_slp *.o
