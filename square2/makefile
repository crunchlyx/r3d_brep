
output: square2.o
	gcc -o square2 square2.o -L../../. -lr3d -lm 

square2.o: square2.c
	gcc -c -g  -std=c99 -I ../../. square2.c

clean:
	rm -r *.o core* square2

