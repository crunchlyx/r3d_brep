
CC = gcc

cube2: cube2.o
	$(CC) -o cube2 cube2.o -L ../../. -lr3d -lm 

cube2.o: cube2.c
	$(CC) -c -g -std=c99 -I ../../. cube2.c

clean:
	rm -r *.o core* cube2
