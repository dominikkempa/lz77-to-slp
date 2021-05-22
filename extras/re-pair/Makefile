SHELL = /bin/sh

CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -DNDEBUG -O3 -pthread -std=c++0x -march=native
CFLAGS_DEBUG = -Wall -Wextra -pedantic -Wshadow -pthread -std=c++0x -g2

all: repair repair_debug

repair:
	$(CC) $(CFLAGS) -o repair ./src/main.cpp ./src/utils.cpp -fopenmp

repair_debug:
	$(CC) $(CFLAGS_DEBUG) -o repair_debug ./src/main.cpp ./src/utils.cpp -fopenmp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f repair repair_debug *.o
