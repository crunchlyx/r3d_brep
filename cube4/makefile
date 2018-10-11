
OPT = -g
CC = gcc
TCMALLOC = -o cube4 cube4.o -L/usr/projects/ngc/private/tcmalloc/gcc-6.4.0/lib -L../../. -lr3d -lm -ltcmalloc
REGMALLOC = -o cube4 cube4.o -L../../. -lr3d -lm

cube4: cube4.o
	$(CC) $(TCMALLOC)

cube4.o: cube4.c
	$(CC) -c $(OPT) -std=c11 -I ../../. cube4.c

clean:
	rm -r *.o core* cube4
