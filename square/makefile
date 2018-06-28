
output: square1.o
	gcc -lm -std=c99 -o square1 square1.o -L../../. -lr3d

square1.o: square1.c
	gcc -c -std=c99 -I ../../. square1.c

clean:
	rm -r *.o core*

