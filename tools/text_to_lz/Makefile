SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -O3 -DNDEBUG -march=native -std=c++0x
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -g2 -std=c++0x

all: text_to_lz

text_to_lz:
	$(CC) $(CFLAGS) -o text_to_lz ./src/main.cpp ./src/utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f text_to_lz *.o
