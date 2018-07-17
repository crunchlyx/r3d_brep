
cube1: cube1.o
	gcc -o cube1 cube1.o -L../../. -lr3d -lm 

cube1.o: cube1.c
	gcc -c -g  -std=c99 -I ../../. cube1.c

clean:
	rm -r *.o core* cube1
