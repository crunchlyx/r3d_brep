
CC = gcc
TCMALLOC = -o cube2 cube2.o -L/usr/projects/ngc/private/tcmalloc/gcc-6.4.0/lib -L../../. -lr3d -lm -ltcmalloc
REGMALLOC = -o cube2 cube2.o -L../../. -lr3d -lm

cube2: cube2.o
	$(CC) $(TCMALLOC)

cube2.o: cube2.c
	$(CC) -c -g -std=c11 -I ../../. cube2.c

clean:
	rm -r *.o core* cube2
