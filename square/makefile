
output: square1.o
	gcc -o square1 square1.o -L../../. -lr3d -lm 

square1.o: square1.c
	gcc -c -O3 -std=c99 -I ../../. square1.c

clean:
	rm -r *.o core* square1

