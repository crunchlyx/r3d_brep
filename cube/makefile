
CC = gcc
TCMALLOC = -o cube1 cube1.o -L/usr/projects/ngc/private/tcmalloc/gcc-6.4.0/lib -L../../. -lr3d -lm -ltcmalloc
REGMALLOC = -o cube1 cube1.o -L../../. -lr3d -lm

cube1: cube1.o
	$(CC) $(TCMALLOC)

cube1.o: cube1.c
	$(CC) -c -g -std=c99 -I ../../. cube1.c

clean:
	rm -r *.o core* cube1
