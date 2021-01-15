SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -O3 -DNDEBUG -march=native
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2

all: test

test:
	$(CC) $(CFLAGS) -o test main.cpp src/utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test *.o
