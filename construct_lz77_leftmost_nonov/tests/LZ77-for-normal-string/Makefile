SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x -march=native
CFLAGS_DEBUG = -Wall -Wextra -Wshadow -std=c++0x -g2

all: test test_debug

test:
	$(CC) $(CFLAGS) -o test main.cpp utils.cpp

test_debug:
	$(CC) $(CFLAGS_DEBUG) -o test_debug main.cpp utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test test_debug *.o
