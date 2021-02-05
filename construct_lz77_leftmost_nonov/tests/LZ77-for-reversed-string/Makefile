SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -Wshadow -funroll-loops -DNDEBUG -O3 -march=native -pthread -std=c++0x
#CFLAGS = -Wall -Wextra -Wshadow -g2 -std=c++0x

all: test

test:
	$(CC) $(CFLAGS) -o test main.cpp utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test *.o
