SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -O3 -DNDEBUG -march=native -std=c++0x
CFLAGS_DEBUG = -Wall -Wextra -pedantic -Wshadow -g2 -std=c++0x

all: test test_debug

test:
	$(CC) $(CFLAGS) -o test main.cpp ../src/utils.cpp

test_debug:
	$(CC) $(CFLAGS_DEBUG) -o test_debug main.cpp ../src/utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test test_debug *.o
