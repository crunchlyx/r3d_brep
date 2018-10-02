
CC = gcc
TCMALLOC = -o cube3 cube3.o -L/usr/projects/ngc/private/tcmalloc/gcc-6.4.0/lib -L../../. -lr3d -lm -ltcmalloc
REGMALLOC = -o cube3 cube3.o -L../../. -lr3d -lm

cube3: cube3.o
	$(CC) $(TCMALLOC)

cube3.o: cube3.c
	$(CC) -c -g -std=c11 -I ../../. cube3.c

clean:
	rm -r *.o core* cube3
