SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -Wshadow -funroll-loops -DNDEBUG -O3 -march=native -pthread -std=c++0x
#CFLAGS = -Wall -Wextra -Wshadow -g2 -std=c++0x

all: bench

bench:
	$(CC) $(CFLAGS) -o bench main.cpp utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f bench *.o
